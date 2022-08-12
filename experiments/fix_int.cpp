#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
int random(int m)
{
  return rand() % m;
}
int main(int argc, const char* argv[])
{
    using namespace Codecset;
    std::string method = std::string(argv[1]);
    std::string source_file = std::string(argv[2]);
    int blocks = atoi(argv[3]);
    int model_size = atoi(argv[4]);

    // We pick a CODEC
    IntegerCODEC &codec = *CODECFactory::getFromName(method);

    std::vector<uint32_t> data;
    std::ifstream srcFile("/home/lyh/Learn-to-Compress/integer_data/"+source_file, std::ios::in);
    // std::ifstream srcFile("/home/lyh/postingList_10M.txt", std::ios::in);
    if (!srcFile)
    {
        
        std::cout << "error opening source file." << std::endl;
        return 0;
    }
    while (srcFile.good())
    {

        uint32_t next;
        srcFile >> next;
        // std::cout<<next<<std::endl;
        if (!srcFile.good())
        {
            break;
        }
        data.push_back(next);
    }
    srcFile.close();
    int N = data.size();
    if (data.size() == 0)
    {
        std::cout << "Empty vector" << std::endl;
        return 0;
    }
    // std::cout << "vector size = " << data.size() << std::endl;
    // std::cout << "vector size = " << data.size() * sizeof(uint32_t) / 1024.0 << "KB"
    //           << std::endl;

    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } 
    // handle with the last block, maybe < block_size
    // std::cout << "Total blocks " << blocks << " block size " << block_size << std::endl;
    int delta = 32;
    codec.init(blocks, block_size, delta);
    std::vector<uint8_t *> block_start_vec;
    std::vector<int> start_index;
    int totalsize = 0;
    double start = getNow();
    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t));
        uint8_t *res = descriptor;
        // std::cout<<data[i*block_size]<<" "<<data[i * block_size+block_size-1]<<std::endl;
        res = codec.encodeArray8(data.data() + (i * block_size), block_length, descriptor, i);
        descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
        block_start_vec.push_back(descriptor);
        totalsize += (res - descriptor);
    }
    double end = getNow();
    double compress_time = end - start;
    double compress_throughput = N*sizeof(uint32_t) / (compress_time*1000000000);

    double origin_size = (sizeof(uint32_t) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size)*100.0 / origin_size;
    double cr_model = total_model_size*100.0 / origin_size;
    double compressrate = (totalsize)*100.0 / origin_size;

    bool flag = true;
    std::vector<uint32_t> recover(data.size());
    double totaltime = 0.0; 
    // std::cout << "decompress all!" << std::endl;
    start = getNow();
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
    end = getNow();
    totaltime += (end - start);
    double da_ns = totaltime / data.size() * 1000000000;

    // std::cout << "all decoding time per int: " << std::setprecision(8)
    //           << totaltime / data.size() * 1000000000 << " ns" << std::endl;
    // std::cout << "all decoding speed: " << std::setprecision(10)
    //           << data.size() / (totaltime * 1000) << std::endl;

    // std::cout << "random access decompress!" << std::endl;
    std::vector<uint32_t> ra_pos;
    int repeat = 1;
    for(int i=0;i<N*repeat;i++){
        ra_pos.push_back(random(N));
    }
    std::vector<uint32_t> buffer(data.size());
    double randomaccesstime = 0.0;
    start = getNow();
    uint32_t mark = 0;
    for (auto index : ra_pos)
    {

        uint32_t tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)index / block_size], index % block_size, buffer.data(), N);
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

    double ra_ns = randomaccesstime / (N*repeat) * 1000000000;

    std::cout<<method<<" "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<cr_model<<" "<<cr_wo_model<<" "<<da_ns<<" "<<ra_ns<<" "<<compress_throughput<<std::endl;

    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
