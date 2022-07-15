// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/leco_string.h"
#include <gperftools/profiler.h>

// typedef leco_uint256 leco_t;
//#include "headers/string/piecewise_fix_string_modify.h"
using namespace Codecset;

template <typename T>
void EncodingOneDataSegment(std::vector<std::string>& data_vec_tmp, int start_ind, int block_length, int padding_length, char padding_char, std::vector<uint8_t*>& descriptor_of_each_block, std::string* common_prefix, int common_prefix_length, uint8_t encoding_type, uint64_t& totalsize_without_padding, uint64_t& totalsize_with_padding, uint64_t& totalsize) {
    Leco_string<T> codec;

    codec.Padding_string(data_vec_tmp, start_ind, block_length,
        padding_char, padding_length);
    codec.get_uncompressed_size(totalsize_with_padding,
        totalsize_without_padding);
    totalsize_with_padding += common_prefix_length * block_length;
    totalsize_without_padding += common_prefix_length * block_length;
    
    uint8_t* descriptor =
        (uint8_t*)malloc(block_length * sizeof(uint64_t) * 1000);
    uint8_t* res = descriptor;

    res = codec.encodeArray8_string(data_vec_tmp, start_ind, block_length,
        descriptor, common_prefix, common_prefix_length, encoding_type);
    int descriptor_length = res - descriptor;
    descriptor = (uint8_t*)realloc(descriptor, descriptor_length);
    totalsize += descriptor_length;
    descriptor_of_each_block.emplace_back(descriptor);
    
    
}


int main()
{

    // ProfilerStart("test_capture.prof");
    std::vector<std::string> string_vec;
    char padding_char = 0;

    // std::ifstream srcFile("/home/lyh/string_data/email_list/mail_server_host_reverse_min_10_max_26.txt", std::ios::in);
    // std::ifstream srcFile("/home/lyh/rocksdb/dump_data/key_2000000_no.txt", std::ios::in);
    std::ifstream srcFile("/home/lyh/Learn-to-Compress/data/string_linear_1M.txt", std::ios::in);
    if (!srcFile)
    {
        std::cout << "error opening source file." << std::endl;
        return 0;
    }
    int cnt = 0;
    while (srcFile.good())
    {
        std::string tmp_str;
        srcFile >> tmp_str;
        if (!srcFile.good())
        {
            break;
        }
        // std::cout << next << std::endl;
        string_vec.push_back(tmp_str);
    }

    srcFile.close();

    std::vector<std::string> string_vec_base(string_vec);
    int N = string_vec.size();
    std::cout << "vector size = " << string_vec.size() << std::endl;

    int block_size = 32;
    int blocks = N / block_size;
    while (block_size * blocks < N)
    {
        blocks++;
    }

    uint64_t totalsize = 0;
    std::vector<uint8_t*> descriptor_of_each_block;
    uint64_t totalsize_without_padding = 0;
    uint64_t totalsize_with_padding = 0;

    for (int i = 0; i < blocks; i++) {
        int block_length = block_size;
        if (i == blocks - 1) {
            block_length = N - (blocks - 1) * block_size;
        }
        int common_prefix_length = 0;
        std::string common_prefix = "";
        int padding_length = extract_common_prefix(
            string_vec, i * block_size, block_length, common_prefix_length, common_prefix);
        // std::cout<<"using common prefix: "<<common_prefix<<" prefix length "<<common_prefix_length<<" padding length: "<<padding_length<<std::endl;
        if (padding_length < 7) {  // the rest can be represented by a uint64_t
            EncodingOneDataSegment<uint64_t>(string_vec, i * block_size, block_length,
                padding_length, padding_char,
                descriptor_of_each_block,
                 &common_prefix, common_prefix_length, 0, totalsize_without_padding, totalsize_with_padding, totalsize);

        }
        else if (padding_length <
            15) {  // the rest can be represented by a uint128_t
            EncodingOneDataSegment<uint128_t>(string_vec, i * block_size, block_length,
                padding_length, padding_char,
                descriptor_of_each_block,
                 &common_prefix, common_prefix_length, 1, totalsize_without_padding, totalsize_with_padding, totalsize);

        }
        else if (padding_length <
            31) {  // the rest can be represented by a uint256_t
            EncodingOneDataSegment<uint256_t>(string_vec, i * block_size, block_length,
                padding_length, padding_char,
                descriptor_of_each_block,
                 &common_prefix, common_prefix_length, 2, totalsize_without_padding, totalsize_with_padding, totalsize);
        }
        else {
            throw("max_padding_length is too large, fall back to mpz");
        }

        // std::cout<<"block "<<i<<" length "<<start_byte<<std::endl;
    }
    totalsize_without_padding+=(N+1);
    totalsize_with_padding+=(N+1);


    double pad_compressrate = (totalsize) * 100.0 / (totalsize_with_padding * 1.0);
    double no_pad_compressrate = (totalsize) * 100.0 / (totalsize_without_padding * 1.0);

    std::cout << N << std::endl;
    std::cout << "compressed size " << totalsize << " without padding uncompressed size " << totalsize_with_padding << " without padding cr%: " << std::setprecision(4) << pad_compressrate << std::endl;
    std::cout << "compressed size " << totalsize << " with padding uncompressed size " << totalsize_without_padding << " with padding cr%: " << std::setprecision(4) << no_pad_compressrate << std::endl;


    bool flag = true;
    double totaltime = 0.0;
    std::cout << "random access decompress!" << std::endl;
    double randomaccesstime = 0.0;
    double start = getNow();
    for (int i = 0; i < N; i++)
    {
        // int index = i;
        int index = random(N);
        std::string result = randomdecode_string(descriptor_of_each_block[index/block_size], index%block_size );
        if (strcmp(result.c_str(), string_vec_base[index].c_str()) != 0)
        {
            flag = false;
            std::cout << "error at index " << index << " result " << result << " expected " << string_vec_base[index] << std::endl;
        }
        if(!flag){
            break;
        }

        
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
