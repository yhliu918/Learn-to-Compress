#include "common.h"
#include "codecfactory.h"
#include "caltime.h"
#include "lr.h"
#include "piecewise_fix_integer_template.h"
#include "piecewise_fix_integer_template_float.h"
#include "piecewise_cost_merge_integer_template.h"
#include "piecewise_cost_merge_integer_template_double.h"
#include "piecewise_cost_merge_hc_integer_template.h"
#include "FOR_integer_template.h"
#include "delta_integer_template.h"
#include "delta_cost_integer_template.h"
#include "delta_cost_merge_integer_template.h"
#include "piecewise_cost_merge_integer_template_test.h"

typedef uint32_t leco_type;

int random(int m)
{
    return rand() % m;
}

template <typename T>
static std::vector<T> load_data_binary(const std::string& filename,
    bool print = true) {
    std::vector<T> data;

    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    // Read size.
    uint64_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    data.resize(size);
    // Read values.
    in.read(reinterpret_cast<char*>(data.data()), size * sizeof(T));
    in.close();
    // std::ofstream out("/home/lyh/Learn-to-Compress/integer_data/wiki_200M_uint64.txt",std::ios::out);
    // std::ios::sync_with_stdio(false);
    // for (auto i = 0; i < size; i++) {
    //   out << data[i] << std::endl;
    // }


    return data;
}


template <typename T>
static std::vector<T> load_data(const std::string& filename) {
    std::vector<T> data;
    std::ifstream srcFile(filename, std::ios::in);
    if (!srcFile) {
        std::cout << "error opening source file." << std::endl;
        return data;
    }

    while (srcFile.good()) {
        T next;
        srcFile >> next;
        if (!srcFile.good()) { break; }
        data.emplace_back(next);

    }
    srcFile.close();

    return data;
}

int main(int argc, const char* argv[])
{
    using namespace Codecset;
    Leco_int<leco_type> codec;
    std::string method = "Leco_fix";
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int delta = atoi(argv[3]);
    int model_size = atoi(argv[4]);
    // alternatives : Delta_int, Delta_cost, Delta_cost_merge, FOR_int, Leco_int, Leco_cost, Leco_cost_merge_hc,  Leco_cost_merge, Leco_cost_merge_double

    // std::vector<leco_type> data = load_data<leco_type>("/home/lyh/Learn-to-Compress/integer_data/" + source_file);
    std::vector<leco_type> data;
    for(int i=0;i<1024*1024*15;i++)
    {
        data.push_back(i);
    }

    int N = data.size();

    // std::cout << "vector size = " << data.size() << std::endl;
    // std::cout << "vector size = " << data.size() * sizeof(leco_type) / 1024.0 << "KB"
    //     << std::endl;

    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size

    // if using auto segmentation codecs
    // int delta = 32;
    // codec.init(blocks, block_size, delta);


    std::vector<uint8_t*> block_start_vec;

    uint64_t totalsize = 0;
    for (int i = 0; i < blocks; i++)
    {
        // std::cout<<"block "<<i<<std::endl;
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        uint8_t* descriptor = (uint8_t*)malloc(block_length * sizeof(leco_type) * 4);
        uint8_t* res = descriptor;
        // if adaptive segment
        // res = codec.encodeArray8_int(data.data(), block_length, descriptor, i);
        // if fixed length segment
        res = codec.encodeArray8_int(data.data() + (i * block_size), block_length, descriptor, i);
        uint32_t segment_size = res - descriptor;
        descriptor = (uint8_t*)realloc(descriptor, segment_size);
        block_start_vec.push_back(descriptor);
        totalsize += segment_size;
    }


    double origin_size = (sizeof(leco_type) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize) * 100.0 / origin_size;

    bool flag = true;
    std::vector<leco_type> recover(data.size());
    double totaltime = 0.0;
    // std::cout << "decompress all!" << std::endl;
    double start = getNow();

    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }

        codec.decodeArray8(block_start_vec[i], block_length, recover.data() + i * block_size, i);

        for (int j = 0; j < block_length; j++)
        {
            if (data[j + i * block_size] != recover[j + i * block_size])
            {
                std::cout << "block: " << i << " num: " << j << " true is: " << data[j + i * block_size] << " predict is: " << recover[j + i * block_size] << std::endl;
                std::cout << "something wrong! decompress failed" << std::endl;
                flag = false;
                break;
            }
        }
        if (!flag)
        {
            break;
        }
    }
    double end = getNow();
    totaltime += (end - start);
    double da_ns = totaltime / data.size() * 1000000000;

    // std::cout << "random access decompress!" << std::endl;
    std::vector<uint32_t> ra_pos;
    int repeat = 1;
    for (int i = 0;i < N * repeat;i++) {
        ra_pos.push_back(random(N));
    }
    flag = true;
    double randomaccesstime = 0.0;
    start = getNow();
    leco_type mark = 0;
    for (auto index : ra_pos)
    {

        leco_type tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)index / block_size], index % block_size, NULL, N);
        // leco_type tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)index / block_size], index, NULL, N);

        mark += tmpvalue;

        if (data[index] != tmpvalue)
        {

            std::cout << "num: " << index << "true is: " << data[index] << " predict is: " << tmpvalue << std::endl;
            flag = false;
            std::cout << "something wrong! decompress failed" << std::endl;
        }
        if (!flag)
        {
            break;
        }
    }
    end = getNow();
    randomaccesstime += (end - start);

    double ra_ns = randomaccesstime / (N * repeat) * 1000000000;


    std::cout << method << " " << source_file << " " << blocks << " " << compressrate << " " << cr_model << " " << cr_wo_model << " " << da_ns << " " << ra_ns << std::endl;


    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
