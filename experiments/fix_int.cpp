#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
int random(int m)
{
    return rand() % m;
}
int main(int argc, const char *argv[])
{
    using namespace Codecset;
    std::string method = std::string(argv[1]);
    std::string source_file = std::string(argv[2]);
    int blocks = atoi(argv[3]);
    int model_size = atoi(argv[4]);
    std::ios::sync_with_stdio(false);
    // We pick a CODEC
    IntegerCODEC &codec = *CODECFactory::getFromName(method);

    std::vector<uint32_t> data;
    std::ifstream srcFile("../data/" + source_file, std::ios::in);

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
    
    // std::sort(data.begin(),data.end());
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    }
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
        uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t)+1000);
        uint8_t *res = descriptor;

        // std::cout<<data[i*block_size]<<" "<<data[i * block_size+block_size-1]<<std::endl;
        res = codec.encodeArray8(data.data() + (i * block_size), block_length, descriptor, i);
        descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
        block_start_vec.push_back(descriptor);
        totalsize += (res - descriptor);
    }
    double end = getNow();
    double compress_time = end - start;

    double compress_throughput = N * sizeof(uint32_t) / (compress_time * 1000000000);

    double origin_size = (sizeof(uint32_t) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize)*100.0 / origin_size;

    bool flag = true;
    std::vector<uint32_t> recover(data.size());
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
            recover[index] += 1;
        }
        for (auto index : codec.mul_add_diff_set_minus)
        {
            recover[index] -= 1;
        }
        for (int j = 0; j < N; j++)
        {
            if (data[j ] != recover[j ])
            {
                std::cout<< " num: " << j << " true is: " << data[j] << " predict is: " << recover[j] << std::endl;
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
    double da_ns = totaltime / (repeat * data.size()) * 1000000000;
    // std::cout << da_ns << std::endl;
    // std::cout << "random access decompress!" << std::endl;
    
    std::vector<uint32_t> ra_pos;
    repeat = 1;
    for(int i=0;i<N*repeat;i++){
        ra_pos.emplace_back(random(N));
        // ra_pos.push_back(i);
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
    std::ofstream outfile("fix_log", std::ios::app);
    outfile<<mark<<std::endl;

    double ra_ns = randomaccesstime / (N*repeat) * 1000000000;

    std::cout<<method<<" "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<cr_model<<" "<<cr_wo_model<<" "<<da_ns<<" "<<ra_ns<<" "<<compress_throughput<<std::endl;
    

    /*
        uint64_t total_sum = 0;
        start = getNow();
        for (int i = 0; i < blocks; i++)
        {
            int block_length = block_size;
            if (i == blocks - 1)
            {
                block_length = N - (blocks - 1) * block_size;
            }
            uint64_t block_sum = codec.summation(block_start_vec[i], block_length, i);
            total_sum += block_sum;
            uint64_t ground_truth = 0;
            for (int j = 0; j < block_length; j++)
            {
                ground_truth+=data[i*block_size+j];
            }
            // std::cout<<"block "<<i<<" "<<ground_truth<<" "<<block_sum<<std::endl;
            if(ground_truth!=block_sum){
                std::cout<<"block "<<i<<" "<<ground_truth<<" "<<block_sum<<std::endl;
            }
            assert(block_sum == ground_truth);
        }
        end = getNow();
        std::cout<<total_sum<<std::endl;
        double sum_time = end - start;
        std::cout << "sum time: " << sum_time/N * 1000000000 << std::endl;
    */
    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
