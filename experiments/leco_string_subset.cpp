// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/leco_string_subset.h"

// typedef leco_uint256 leco_t;
//#include "headers/string/piecewise_fix_string_modify.h"
using namespace Codecset;

// email: 46-123
//uuid: 46-90
// words: 97 - 122
int max_ascii = 122;
int min_ascii = 97;
template <typename T>
void EncodingOneDataSegment(std::vector<std::string>& data_vec_tmp, int start_ind, int block_length, int padding_length, char padding_char, std::vector<uint8_t*>& descriptor_of_each_block, std::string* common_prefix, int common_prefix_length, uint8_t encoding_type, uint64_t& totalsize_without_padding, uint64_t& totalsize_with_padding, uint64_t& totalsize) {
    Leco_string_subset<T> codec;
    codec.set_ascii_set(min_ascii, max_ascii);
    codec.Padding_string(data_vec_tmp, start_ind, block_length,
        padding_char, padding_length, min_ascii, max_ascii);
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




int main(int argc, const char* argv[])
{
    std::string source_file = std::string(argv[1]);
    int block_size = atoi(argv[2]);

    std::vector<std::string> string_vec;
    char padding_char = min_ascii;
    std::ifstream srcFile(source_file, std::ios::in);
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
    // std::string tmp_str;
    // for(std::string line;std::getline(srcFile, tmp_str, srcFile.widen('\n'));){
    //     string_vec.push_back(tmp_str);
    //     uint128_t result = convertToASCII<uint128_t>(tmp_str);
    //     print_u128_u(result);
    //     std::cout<<std::endl;
    // }

    

    std::vector<std::string> string_vec_base(string_vec);
    int N = string_vec.size();
    std::cout << "vector size = " << string_vec.size() << std::endl;

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
    uint64_t totalsize_with_index = totalsize+blocks*sizeof(unsigned);
    double compressrate_with_index = totalsize_with_index*100.0 / (totalsize_without_padding * 1.0);

    std::cout << N << std::endl;
    std::cout << "compressed size " << totalsize << " with padding uncompressed size " << totalsize_without_padding << " with padding cr%: " << std::setprecision(4) << no_pad_compressrate << std::endl;
    std::cout << "compressed size with index " << totalsize_with_index << " with padding uncompressed size " << totalsize_without_padding << " with padding cr%: " << std::setprecision(4) << compressrate_with_index << std::endl;

    std::vector<int> ra_pos;
    int repeat = 100;
    for(int i=0;i<N*repeat;i++){
        ra_pos.push_back(random(N));
        // ra_pos.push_back(i);
    }

    bool flag = true;
    double totaltime = 0.0;
    std::cout << "random access decompress!" << std::endl;
    double randomaccesstime = 0.0;
    double start = getNow();
    
    for (auto index : ra_pos)
    {

        std::string result = randomdecode_string(descriptor_of_each_block[index/block_size], index%block_size,min_ascii, max_ascii );
        // std::string result = randomdecode_string(descriptor_of_each_block[index/block_size], index%block_size);
        // if (memcmp(result.c_str(), string_vec_base[index].c_str(), string_vec_base[index].size()) != 0)
        // {
        //     flag = false;
        //     std::cout << "error at index " << index << " result " << result << " expected " << string_vec_base[index] << std::endl;
        // }
        // if(!flag){
        //     break;
        // }

        
    }
    double end = getNow();
    randomaccesstime += (end - start)/repeat;

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
