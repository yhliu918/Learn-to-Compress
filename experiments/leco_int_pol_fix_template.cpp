#include "common.h"
#include "codecfactory.h"
#include "caltime.h"
#include "lr.h"
#include "poly_fix_integer_template.h"
#include "poly_fix_integer_template_max.h"
#include "piecewise_fix_integer_template_float.h"
#include "piecewise_cost_merge_integer_template_double.h"
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
    
    std::string method = "Poly_fix_op_max";
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int model_size = atoi(argv[3]);
    bool binary = atoi(argv[4]);
    leco_type filter1 = 0;
    leco_type filter2 = 0;
    leco_type base = 0;
    bool filter_experiment = false;
    bool filter_close_experiment = false;
    if (argc > 5)
    {
        filter1 = atoll(argv[5]);
        filter_experiment = true;
    }
    if (argc > 6)
    {
        filter2 = atoll(argv[6]);
        filter_experiment = false;
        filter_close_experiment = true;
        base = atoll(argv[7]);
    }
    // alternatives : Delta_int, Delta_cost, Delta_cost_merge, FOR_int, Leco_int, Leco_cost, Leco_cost_merge_hc,  Leco_cost_merge, Leco_cost_merge_double

    std::vector<leco_type> data;
    if(binary){
        data = load_data_binary<leco_type>(source_file);
        // data = load_data_binary<leco_type>("../data/" + source_file);
    }
    else{
        data = load_data<leco_type>(source_file);
        // data = load_data<leco_type>("../data/" + source_file);
    }
    int N = data.size();

    // std::string seg_file = "/root/Learn-to-Compress/buildnew/pol_seg.log";
    // std::vector<int> seg_idx = load_data<int>(seg_file);
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size
    // convert block_size to const
    // const int contsblock_size = const_cast<const int&>(block_size);
    Leco_int_poly<leco_type, 3> codec;
    codec.init(blocks, block_size);

    std::vector<uint8_t*> block_start_vec;
    double start = getNow();
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
        res = codec.encodeArray8_int(data.data() + (i * block_size), block_length, descriptor, i);
        uint32_t segment_size = res - descriptor;
        descriptor = (uint8_t*)realloc(descriptor, segment_size);
        block_start_vec.push_back(descriptor);
        totalsize += segment_size;
    }
    // for (int i = 0; i < seg_idx.size()-1; i++)
    // {
    //     // std::cout<<"block "<<i<<std::endl;
    //     int block_length = seg_idx[i+1] - seg_idx[i];

    //     uint8_t* descriptor = (uint8_t*)malloc(block_length * sizeof(leco_type) * 4+1000);
    //     uint8_t* res = descriptor;
    //     res = codec.encodeArray8_int(data.data() + seg_idx[i], block_length, descriptor, i);
    //     uint32_t segment_size = res - descriptor;
    //     descriptor = (uint8_t*)realloc(descriptor, segment_size);
    //     block_start_vec.push_back(descriptor);
    //     totalsize += segment_size;
    // }
    double end = getNow();
    double compress_time = end - start;
    double compress_throughput = N * sizeof(leco_type) / (compress_time * 1000000000);


    double origin_size = (sizeof(leco_type) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize) * 100.0 / origin_size;

    bool flag = true;
    std::vector<leco_type> recover(data.size());
    double totaltime = 0.0;
    // std::cout << "decompress all!" << std::endl;
    int repeat = 10;

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
            codec.decodeArray8(block_start_vec[i], block_length, recover.data() + i * block_size, i);
        }
        for (auto index : codec.mul_add_diff_set)
        {
            recover[index.first] += index.second;
        }
        // #ifndef NDEBUG
        // for (int j = 0; j < N; j++)
        // {
        //     if (data[j ] != recover[j ])
        //     {
        //         std::cout<< " num: " << j << " true is: " << data[j] << " predict is: " << recover[j] << std::endl;
        //         std::cout << "something wrong! decompress all failed" << std::endl;
        //         flag = false;
        //         break;
        //     }
        // }
        // if (!flag)
        // {
        //     break;
        // }
        // #endif
    }
    end = getNow();
    totaltime += (end - start);
    double da_ns = totaltime / (N*repeat) * 1000000000;

    // std::cout << "random access decompress!" << std::endl;
    std::vector<uint32_t> ra_pos;
    repeat = 1;
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
        mark += tmpvalue;
        // #ifndef NDEBUG
        // if (data[index] != tmpvalue)
        // {

        //     std::cout << "num: " << index << "true is: " << data[index] << " predict is: " << tmpvalue << std::endl;
        //     flag = false;
        //     std::cout << "something wrong! random access failed" << std::endl;
        // }
        // if (!flag)
        // {
        //     break;
        // }
        // #endif
    }
    end = getNow();
    randomaccesstime += (end - start);
    std::ofstream outfile("fix_log", std::ios::app);
    outfile<<mark<<std::endl;

    double ra_ns = randomaccesstime / (N * repeat) * 1000000000;


    std::cout << method << " " << source_file << " " << blocks << " " << compressrate << " " << cr_model << " " << cr_wo_model << " " << da_ns << " " << ra_ns<<" "<<compress_throughput<< std::endl;


    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
