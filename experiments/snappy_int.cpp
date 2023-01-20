#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
#include "snappy.h"
typedef uint64_t leco_type;
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

int main(int argc, const char *argv[])
{
    using namespace Codecset;
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int model_size = atoi(argv[3]);
    bool binary = atoi(argv[4]);
    std::ios::sync_with_stdio(false);
    // We pick a CODEC

    std::vector<leco_type> data;
    source_file = "../data/"+source_file;

    if(binary){
        data = load_data_binary<leco_type>(source_file);
    }
    else{
        data = load_data<leco_type>(source_file);
    }
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
    std::vector<std::string> block_start_vec;
    std::vector<int> start_index;
    std::vector<int> seglen;

    int totalsize = 0;
    double start = getNow();
    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        // uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t)+1000);
        // uint8_t *res = descriptor;
        std::string compressed;
        int segment_size = snappy::Compress((char *)(data.data() + i * block_size), block_length * sizeof(leco_type), &compressed);
        // std::cout<<segment_size<<std::endl;
        seglen.push_back(segment_size);
        // descriptor = (uint8_t *)realloc(descriptor,segment_size);
        block_start_vec.push_back(compressed);
        totalsize += segment_size;
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
            std::string uncomp_str;
            // const char* reover= reinterpret_cast<const char*>(recover + block_size*i);
            snappy::Uncompress(block_start_vec[i].c_str(), seglen[i], &uncomp_str);
            // std::cout<<uncomp_str<<std::endl;
            // &recover[i*block_size] = reinterpret_cast<leco_type*>
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
        std::string uncomp_str;
        snappy::Uncompress(block_start_vec[(int)index / block_size].c_str(), seglen[(int)index / block_size], &uncomp_str);
        leco_type tmpvalue;
        memcpy(&tmpvalue, uncomp_str.c_str() + (index % block_size) * sizeof(leco_type), sizeof(leco_type));
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
    outfile<<mark<<std::endl;

    double ra_ns = randomaccesstime / (N*repeat) * 1000000000;

    std::cout<<"snappy "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<cr_model<<" "<<cr_wo_model<<" "<<da_ns<<" "<<ra_ns<<" "<<compress_throughput<<std::endl;
    

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
    // for (int i = 0; i < (int)block_start_vec.size(); i++)
    // {
    //     free(block_start_vec[i]);
    // }
}
