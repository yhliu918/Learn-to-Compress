/**
 * This is an implementation of Huffman coding.
 *
 * The core algorithm is taken from the CLR book (Introduction of Algorithms),
 * Chapter 16.3, and directly used to implement the 'build_tree()' routine.
 *
 * After the tree is built, a code table that maps a character to a binary
 * code is built from the tree, and used for encoding text. Decoding is done
 * by traversing the Huffman tree, as prescribed by the algorithm.
 *
 * Binary codes are represented by std::vector<bool>, which is a specialized
 * vector that optimizes space.
 */

#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <string>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include "../headers/common.h"
#include "../headers/platform.h"
#include "../headers/rans_byte.h"
typedef uint32_t leco_type;
using namespace std;

// ---- Stats

struct SymbolStats
{
    uint32_t freqs[256];
    uint32_t cum_freqs[257];

    void count_freqs(uint8_t const *in, size_t nbytes);
    void calc_cum_freqs();
    void normalize_freqs(uint32_t target_total);
};

void SymbolStats::count_freqs(uint8_t const *in, size_t nbytes)
{
    for (int i = 0; i < 256; i++)
        freqs[i] = 0;

    for (size_t i = 0; i < nbytes; i++)
        freqs[in[i]]++;
}

void SymbolStats::calc_cum_freqs()
{
    cum_freqs[0] = 0;
    for (int i = 0; i < 256; i++)
        cum_freqs[i + 1] = cum_freqs[i] + freqs[i];
}

void SymbolStats::normalize_freqs(uint32_t target_total)
{
    assert(target_total >= 256);

    calc_cum_freqs();
    uint32_t cur_total = cum_freqs[256];

    // resample distribution based on cumulative freqs
    for (int i = 1; i <= 256; i++)
        cum_freqs[i] = ((uint64_t)target_total * cum_freqs[i]) / cur_total;

    // if we nuked any non-0 frequency symbol to 0, we need to steal
    // the range to make the frequency nonzero from elsewhere.
    //
    // this is not at all optimal, i'm just doing the first thing that comes to mind.
    for (int i = 0; i < 256; i++)
    {
        if (freqs[i] && cum_freqs[i + 1] == cum_freqs[i])
        {
            // symbol i was set to zero freq

            // find best symbol to steal frequency from (try to steal from low-freq ones)
            uint32_t best_freq = ~0u;
            int best_steal = -1;
            for (int j = 0; j < 256; j++)
            {
                uint32_t freq = cum_freqs[j + 1] - cum_freqs[j];
                if (freq > 1 && freq < best_freq)
                {
                    best_freq = freq;
                    best_steal = j;
                }
            }
            assert(best_steal != -1);

            // and steal from it!
            if (best_steal < i)
            {
                for (int j = best_steal + 1; j <= i; j++)
                    cum_freqs[j]--;
            }
            else
            {
                assert(best_steal > i);
                for (int j = i + 1; j <= best_steal; j++)
                    cum_freqs[j]++;
            }
        }
    }

    // calculate updated freqs and make sure we didn't screw anything up
    assert(cum_freqs[0] == 0 && cum_freqs[256] == target_total);
    for (int i = 0; i < 256; i++)
    {
        if (freqs[i] == 0)
            assert(cum_freqs[i + 1] == cum_freqs[i]);
        else
            assert(cum_freqs[i + 1] > cum_freqs[i]);

        // calc updated freq
        freqs[i] = cum_freqs[i + 1] - cum_freqs[i];
    }
}

int random(int m)
{
    return rand() % m;
}
template <typename T>
static std::vector<T> load_data_binary(const std::string &filename,
                                       bool print = true)
{
    std::vector<T> data;

    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open())
    {
        std::cerr << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    // Read size.
    uint64_t size;
    in.read(reinterpret_cast<char *>(&size), sizeof(uint64_t));
    data.resize(size);
    // Read values.
    in.read(reinterpret_cast<char *>(data.data()), size * sizeof(T));
    in.close();

    return data;
}

template <typename T>
static std::vector<T> load_data(const std::string &filename)
{
    std::vector<T> data;
    std::ifstream srcFile(filename, std::ios::in);
    if (!srcFile)
    {
        std::cout << "error opening source file." << std::endl;
        return data;
    }

    while (srcFile.good())
    {
        T next;
        srcFile >> next;
        if (!srcFile.good())
        {
            break;
        }
        data.emplace_back(next);
    }
    srcFile.close();

    return data;
}
double getNow()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}
static const uint32_t prob_bits = 14;
static const uint32_t prob_scale = 1 << prob_bits;

int main(int argc, const char *argv[])
{
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int model_size = atoi(argv[3]);
    bool binary = atoi(argv[4]);
    std::vector<leco_type> data;
    source_file = "../data/" + source_file;

    if (binary)
    {
        data = load_data_binary<leco_type>(source_file);
    }
    else
    {
        data = load_data<leco_type>(source_file);
    }
    int N = data.size();
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    }
    int delta = 32;
    std::vector<uint8_t*> block_start_vec;
    std::vector<int> start_index;
    std::vector<int> seglen;
    std::vector<uint8_t *> rans_begins;
    std::vector<uint8_t*> cumsums;
    std::vector<RansDecSymbol*> decsym;

    int totalsize = 0;
    double start = getNow();
    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        uint8_t *in_bytes = reinterpret_cast<uint8_t *>((data.data() + i * block_size));
        int in_size = block_length * sizeof(leco_type);
        SymbolStats stats;
        stats.count_freqs(in_bytes, in_size);
        stats.normalize_freqs(prob_scale);

        uint8_t* cum2sym = new uint8_t[prob_scale];
        for (int s = 0; s < 256; s++)
            for (uint32_t i = stats.cum_freqs[s]; i < stats.cum_freqs[s + 1]; i++)
                cum2sym[i] = s;
        cumsums.push_back(cum2sym);
        static size_t out_max_size = 32 << 20; // 32MB
        uint8_t *out_buf = new uint8_t[out_max_size];

        uint8_t *rans_begin;
        RansEncSymbol esyms[256];
        RansDecSymbol* dsyms = new RansDecSymbol[256];

        for (int i = 0; i < 256; i++)
        {
            RansEncSymbolInit(&esyms[i], stats.cum_freqs[i], stats.freqs[i], prob_bits);
            RansDecSymbolInit(&dsyms[i], stats.cum_freqs[i], stats.freqs[i]);
        }

        decsym.push_back(dsyms);

        RansState rans;
        RansEncInit(&rans);

        uint8_t *ptr = out_buf + out_max_size; // *end* of output buffer
        for (size_t i = in_size; i > 0; i--)
        { // NB: working in reverse!
            int s = in_bytes[i - 1];
            RansEncPutSymbol(&rans, &ptr, &esyms[s]);
        }
        RansEncFlush(&rans, &ptr);
        rans_begin = ptr;

        // printf("rANS: %d bytes\n", (int)(out_buf + out_max_size - rans_begin));
        rans_begins.push_back(rans_begin);
        block_start_vec.push_back(out_buf);
        totalsize += (out_buf + out_max_size - rans_begin +  256*4 );
    }
    double end = getNow();
    double compress_time = end - start;

    double compress_throughput = N * sizeof(leco_type) / (compress_time * 1000000000);

    double origin_size = (sizeof(leco_type) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize)*100.0 / origin_size;

    bool flag = true;
    std::vector<leco_type> recover(data.size());
    double totaltime = 0.0;
    // std::cout << "decompress all!" << std::endl;
    int repeat = 1;

    start = getNow();
    for (int iter = 0; iter < repeat; iter++)
    {
        for (int i = 0; i < blocks; i++)
        {
            int block_length = block_size;
            if (i == blocks - 1)
            {
                block_length = N - (blocks - 1) * block_size;
            }
            RansState rans;
            uint8_t* dec_bytes = reinterpret_cast<uint8_t*> (recover.data()+i*block_size);
            uint8_t *ptr = rans_begins[i];
            RansDecInit(&rans, &ptr);
            int in_size = block_length*sizeof(leco_type);
            RansDecSymbol* dsyms = decsym[i];
            uint8_t* cumsum = cumsums[i];
            for (size_t i = 0; i < in_size; i++)
            {
                uint32_t s = cumsum[RansDecGet(&rans, prob_bits)];
                dec_bytes[i] = (uint8_t)s;
                RansDecAdvanceSymbol(&rans, &ptr, &dsyms[s], prob_bits);
            }

        }
    }
    end = getNow();

    totaltime += (end - start);
    double da_ns = totaltime / (repeat * data.size()) * 1000000000;
    // std::cout << da_ns << std::endl;
    // std::cout << "random access decompress!" << std::endl;

    std::vector<uint32_t> ra_pos;
    // repeat = 1;
    int total_ra = N/1000;
    for (int i = 0; i < total_ra; i++)
    {
        ra_pos.emplace_back(random(N));
        // ra_pos.push_back(i);
    }
    std::vector<uint32_t> buffer(data.size());

    double randomaccesstime = 0.0;
    uint32_t mark = 0;
    int in_size = block_size*sizeof(leco_type);
    uint8_t* dec_bytes = new uint8_t[block_size*sizeof(leco_type)];
    start = getNow();
    for (auto index : ra_pos)
    {
        RansState rans;
        uint8_t *ptr = rans_begins[index/block_size];
        RansDecInit(&rans, &ptr);
        RansDecSymbol* dsyms = decsym[index/block_size];
        uint8_t* cumsum = cumsums[index/block_size];
        for (size_t i = 0; i < in_size; i++)
        {
            uint32_t s = cumsum[RansDecGet(&rans, prob_bits)];
            dec_bytes[i] = (uint8_t)s;
            RansDecAdvanceSymbol(&rans, &ptr, &dsyms[s], prob_bits);
        }

        leco_type tmpvalue = reinterpret_cast<leco_type*>(dec_bytes)[index%block_size];
        
        // std::cout<<tmpvalue<<std::endl;
        mark += tmpvalue;
        // if (data[index] != tmpvalue)
        // {

        //     std::cout << "num: " << index << "true is: " << data[index] << " predict is: " << tmpvalue << std::endl;
        //     flag = false;
        //     std::cout << "something wrong! decompress failed" << std::endl;
        // }
        // if (!flag)
        // {
        //     break;
        // }
    }

    end = getNow();
    randomaccesstime += (end - start);
    std::ofstream outfile("fix_log", std::ios::app);
    outfile << mark << std::endl;

    double ra_ns = randomaccesstime / total_ra * 1000000000;

    std::cout << "RNS " << source_file << " " << blocks << " " << compressrate << " " << cr_model << " " << cr_wo_model << " " << da_ns << " " << ra_ns << " " << compress_throughput << std::endl;
}
