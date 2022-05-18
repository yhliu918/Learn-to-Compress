// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/piecewise_fix_string_padding_template.h"
#include <gperftools/profiler.h>
#include "../headers/leco_uint256.h"
//typedef uint128_t leco_t;
typedef leco_uint256 leco_t;
//#include "headers/string/piecewise_fix_string_modify.h"
using namespace Codecset;


int main()
{

    // ProfilerStart("test_capture.prof");
    std::vector<std::string> string_vec;
    Leco_string<leco_t> codec;
    char padding_char = 2;

    std::ifstream srcFile("/home/lyh/rocksdb/dump_data/wholestring_max_30_sep_key.txt", std::ios::in);
    // std::ifstream srcFile("/home/lyh/rocksdb/dump_data/poisson_20000_sep_key.txt", std::ios::in);
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

    
    int N = string_vec.size();
    
    int block_size = 64;
    int blocks = N / block_size;
    while (block_size * blocks < N)
    {
        blocks++;
    }
    codec.setBlockSize(block_size, blocks);
    codec.Padding_string(string_vec, N, padding_char);

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
        std::string tmp_str = string_vec[i * block_size];
        codec.Padd_one_string(tmp_str, padding_char);
        codec.push_back_firstkey(tmp_str);
    }
    int max_length = codec.get_max_padding_length();
    int with_padding = 0;
    int without_padding = 0;
    int whole_padding_size  = N * max_length;
    codec.get_uncompressed_size(with_padding, without_padding);
    int overhead = blocks * max_length + N; // first key and string length(each uint8_t)
    totalsize += overhead;

    double pad_compressrate = (totalsize)*100.0 / (with_padding * 1.0);
    double no_pad_compressrate = (totalsize)*100.0 / (without_padding * 1.0);
    double whole_pad_compressrate = (totalsize)*100.0 / (whole_padding_size * 1.0);
    
    std::cout << N <<std::endl;
    std::cout << "compressed size " << totalsize << " without padding uncompressed size " << without_padding <<" without padding cr%: "<<std::setprecision(4)<<no_pad_compressrate<< std::endl;
    std::cout << "compressed size " << totalsize << " with padding uncompressed size " << with_padding <<" with padding cr%: "<<std::setprecision(4)<<pad_compressrate<< std::endl;
    std::cout << "compressed size " << totalsize << " whole padding uncompressed size " << whole_padding_size <<" with padding cr%: "<<std::setprecision(4)<<whole_pad_compressrate<< std::endl;


    bool flag = true;
    double totaltime = 0.0;
    std::cout << "random access decompress!" << std::endl;
    double randomaccesstime = 0.0;
    double start = getNow();
    for (int i = 0; i < N; i++)
    {
        int index = i;
        leco_t tmpval;

        // std::string tmpstr = codec.randomdecodeArray8_string((int)index / block_size, index % block_size, &tmpval);

        // std::cout<<index<<" "<<index % block_size<<" block "<<(int)index / block_size<<std::endl;
        // std::cout << std::endl;

        
        // assert(tmpstr.size() !=0);
        // std::string target = string_vec[index].substr(0, tmpstr.size());
        // assert(tmpstr == target);

        // if (target.compare(tmpstr))
        // {
        // std::cout << "num: " << index<<" block "<< index/block_size << " true is: " << target << " predict is: " << tmpstr << std::endl;
        // flag = false;

        // std::cout << "something wrong! decompress failed" << std::endl;
        // }
        // if (!flag)
        // {
        // break;
        // }
    }
    double end = getNow();
    randomaccesstime += (end - start);

    std::cout << "random decoding time per int: " << std::setprecision(8)
              << randomaccesstime / N * 1000000000 << " ns" << std::endl;


    // int sample_size = 1000;
    // start = getNow();
    // codec.TestBsearch(sample_size, string_vec, N);
    // end = getNow();
    // double ourbinarytime = end - start;
    // std::cout << "binary time per time: " << std::setprecision(8)
    //           << ourbinarytime / sample_size * 1000000000 << " ns" << std::endl;

    // ProfilerStop();
}
