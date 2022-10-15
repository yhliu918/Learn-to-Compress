// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/caltime.h"
#include "../headers/string/fsst_string.h"
#include "../headers/bit_read.h"
#include "../headers/bit_write.h"
#include "../headers/delta_my.h"
using namespace Codecset;

int random(int m)
{
  return rand() % m;
}
int main(int argc, const char* argv[])
{
    std::string source_file = std::string(argv[1]);
    bool offset_delta_compress = std::atoi(argv[2]);
    int delta_block_size = std::atoi(argv[3]);

    std::vector<std::string> string_vec;
    char padding_char = 0;
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

    FSST_string codec;
    delta_my delta_codec;
    

    std::vector<std::string> string_vec_base(string_vec);
    int N = string_vec.size();
    int block_size = N;
    // std::cout << "vector size = " << string_vec.size() << std::endl;

    int blocks = N / block_size;
    while (block_size * blocks < N)
    {
        blocks++;
    }
    int delta_blocks = N/delta_block_size;
    while (delta_block_size * delta_blocks < N)
    {
        delta_blocks++;
    }
    if(offset_delta_compress){
        delta_codec.init(delta_blocks, delta_block_size, 0);
    }

    uint64_t totalsize = 0;
    uint64_t totalsize_delta_offset = 0;
    std::vector<uint8_t*> descriptor_of_each_block;
    uint64_t totalsize_without_padding = 0;
    std::vector<uint32_t> offset;
    std::vector<uint8_t*> delta_descriptor_of_each_block;


    for(auto item: string_vec)
    {
        totalsize_without_padding += item.size();
    }
    for (int i = 0; i < blocks; i++) {
        int block_length = block_size;
        if (i == blocks - 1) {
            block_length = N - (blocks - 1) * block_size;
        }
        totalsize = codec.encodeArray8_string(string_vec, offset);
    }
    offset.insert(offset.begin(), 0);
    if(offset_delta_compress){
        for (int i = 0; i < delta_blocks; i++) {
            int block_length = delta_block_size;
            if (i == delta_blocks - 1) {
                block_length = N - (delta_blocks - 1) * delta_block_size;
            }
            uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t)*2);
            uint8_t *res = descriptor;
            res = delta_codec.encodeArray8(offset.data() + (i * delta_block_size), block_length, descriptor, i);
            descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
            delta_descriptor_of_each_block.push_back(descriptor);
            totalsize_delta_offset += (res - descriptor);
        }
        
    }

    totalsize_without_padding+=(N+1);
    double no_pad_compressrate = (totalsize) * 100.0 / (totalsize_without_padding * 1.0);

    uint64_t totalsize_with_index = N*sizeof(uint32_t)+totalsize;
    if(offset_delta_compress){
        totalsize_with_index = totalsize_delta_offset+totalsize;
    }
    double no_pad_compressrate_with_ind = (totalsize_with_index)*100.0/(totalsize_without_padding*1.0);

    std::vector<int> ra_pos;
    int repeat = 100;
    for(int i=0;i<N*repeat;i++){
        ra_pos.push_back(random(delta_block_size*(N/delta_block_size)));
        // ra_pos.push_back(i);
    }

    bool flag = true;
    double totaltime = 0.0;
    // std::cout << "random access decompress!" << std::endl;
    double randomaccesstime = 0.0;
    double start = getNow();
    std::vector<uint32_t> buffer(N);
    uint32_t* out = new uint32_t[delta_block_size];
    uint32_t* out2 = new uint32_t[delta_block_size];
    for (auto index : ra_pos)
    {
        std::string result;
        if(offset_delta_compress){
            
            uint32_t offset_val = 0;
            uint32_t offset_val2 = 0;

            out = delta_codec.decodeArray8(delta_descriptor_of_each_block[(int)(index) / delta_block_size], delta_block_size,out, N);
            offset_val = out[(index) % delta_block_size];
            
            if((index+1)%delta_block_size==0 && (index+1)<N){
                out2 = delta_codec.decodeArray8(delta_descriptor_of_each_block[(int)(index+1) / delta_block_size], delta_block_size,out2, N);
                offset_val2 = out2[(index+1)%delta_block_size];
            }
            else{
                offset_val2 = out[(index+1)%delta_block_size];
            }
            if(index==N-1){
                offset_val2 = offset_val+100;
            }
            // uint32_t offset_val2 = delta_codec.randomdecodeArray8(delta_descriptor_of_each_block[index/delta_block_size], index%delta_block_size, out, N);
            // std::cout<<index<<" "<<offset_val<<" "<<offset_val2<<std::endl;
            result = codec.randomdecode_string(index, offset_val, offset_val2);
            

        }
        else{
            uint32_t offset_val=0;
            if(index){
                offset_val = offset[(index-1)];
            }
            result = codec.randomdecode_string(index, offset_val, offset[index]);
        }
        
        // if (memcmp(result.c_str(), string_vec[index].c_str(), string_vec[index].size()) != 0)
        // {
        //     flag = false;
        //     std::cout << "error at index " << index << " result " << result << " expected " << string_vec[index] << std::endl;
        // }
        // if(!flag){
        //     break;
        // }

        
    }
    double end = getNow();
    randomaccesstime += (end - start)/repeat;
    double ra_ns = randomaccesstime / N * 1000000000;

    // std::cout << "random decoding time per int: " << std::setprecision(8)
    //     << randomaccesstime / N * 1000000000 << " ns" << std::endl;

    std::cout<<source_file<<" "<<delta_block_size<<" "<<no_pad_compressrate<<" "<<no_pad_compressrate_with_ind<<" "<<ra_ns<<std::endl;

    // int sample_size = 1000;
    // start = getNow();
    // codec.TestBsearch(sample_size, string_vec, N);
    // end = getNow();
    // double ourbinarytime = end - start;
    // std::cout << "binary time per time: " << std::setprecision(8)
    //           << ourbinarytime / sample_size * 1000000000 << " ns" << std::endl;
    free(out);
    free(out2);
    // ProfilerStop();
}
