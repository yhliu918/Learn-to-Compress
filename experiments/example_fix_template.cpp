// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/piecewise_fix_string_template.h"
#include <gperftools/profiler.h>
#include "../headers/leco_uint256.h"
//typedef leco_uint256 leco_t;
typedef uint128_t leco_t;
//#include "headers/string/piecewise_fix_string_modify.h"
using namespace Codecset;


int main()
{

    // ProfilerStart("test_capture.prof");
    std::vector<std::string> string_vec;
    Leco_string<leco_t> codec;

    // std::ifstream srcFile("/home/lyh/rocksdb/dump_data/padding_a_wholestring_max_30.txt", std::ios::in);
    std::ifstream srcFile("/home/lyh/Learn-to-Compress/data/poisson_2000000/key.txt", std::ios::in);
    if (!srcFile)
    {
        std::cout << "error opening source file." << std::endl;
        return 0;
    }
    int cnt = 0;
    while (1)
    {
        std::string tmp_str;
        srcFile >> tmp_str;
        if (srcFile.eof())
        {
            break;
        }
        // std::cout << next << std::endl;
        string_vec.push_back(tmp_str);
    }

    srcFile.close();

    int string_length = string_vec[0].size();
    int N = string_vec.size();
    int blocks = 200;
    int block_size = N / blocks;
    codec.setBlockSize(block_size);
    if (block_size * blocks < N)
    {
        blocks++;
    }

    uint64_t totalsize = 0;
    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }

        uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t) * 1000);
        uint8_t *res = descriptor;

        res = codec.encodeArray8_string(string_vec, i * block_size, block_length, descriptor, i);
        totalsize += (res - descriptor);
        descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
        codec.push_back_block(descriptor);
    }

    for (int i = 0; i < blocks; i++)
    {
        codec.push_back_firstkey(string_vec[i * block_size]);
    }
    totalsize += blocks * 30;
    double compressrate = (totalsize)*100.0 / (N * string_length * 1.0);
    std::cout << N << " " << string_length << " " << totalsize << " " << compressrate << std::endl;
    std::cout << "compressed size " << totalsize << " uncompressed size " << N * string_length << std::endl;
    std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;


    bool flag = true;
    double totaltime = 0.0;
    std::cout << "random access decompress!" << std::endl;
    double randomaccesstime = 0.0;
    double start = getNow();
    for (int i = 0; i < N; i++)
    {
        int index = random(N);
        leco_t tmpval;
        codec.randomdecodeArray8((int)index / block_size, index % block_size, &tmpval);
        //std::string tmp_str = convertToString<leco_t>(&tmpval);
        // // std::cout<<i<<" "<<tmpvalue<<" block "<<(int)index / block_size<<std::endl;
        //  if (string_vec[index].compare(tmp_str))
        //  {
        //    std::cout << "num: " << index<<" block "<< index/block_size << " true is: " << string_vec[index] << " predict is: " << tmp_str << std::endl;
        //    flag = false;
        //    std::cout << "something wrong! decompress failed" << std::endl;
        //  }
        //  if (!flag)
        //  {
        //    break;
        //  }
    }
    double end = getNow();
    randomaccesstime += (end - start);

    std::cout << "random decoding time per int: " << std::setprecision(8)
              << randomaccesstime / N * 1000000000 << " ns" << std::endl;


    int sample_size = 1000;
    start = getNow();
    codec.TestBsearch(sample_size, string_vec, N);
    end = getNow();
    double ourbinarytime = end - start;
    std::cout << "binary time per time: " << std::setprecision(8)
              << ourbinarytime / sample_size * 1000000000 << " ns" << std::endl;

    // ProfilerStop();
}
