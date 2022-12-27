#include <iostream>
#include <vector>
#include <iomanip>

#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "util.hpp"

#include "elias_fano.hpp"
#include "mapper.hpp"

#include "perftest_common.hpp"

typedef uint64_t data_type;

int random(int m)
{
  return rand() % m;
}
double getNow() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1000000.0;
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
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int binary = atoi(argv[3]);

    std::vector<data_type> data;
    if(!binary){
        data = load_data<data_type>("/root/Learn-to-Compress/data/" + source_file);
    }
    else{
        data = load_data_binary<data_type>("/root/Learn-to-Compress/data/" + source_file);
    }

    int N = data.size();
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size
    
    int block_length = block_size;
    
    std::vector<uint8_t*> block_start_vec;
    uint64_t totalsize = 0;
    double start_cr = getNow();
    for(int i=0;i<blocks;i++)
    {
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        uint64_t max_element = data[i*block_size+block_length -1];
        succinct::elias_fano::elias_fano_builder* tmp_bvb_size = new succinct::elias_fano::elias_fano_builder(max_element, block_length);
        for(int j=0;j<block_length;j++)
        {
            (*tmp_bvb_size).push_back(data[i*block_size+j]);
        } 
        succinct::elias_fano* ef = new succinct::elias_fano(tmp_bvb_size);
        uint8_t* descriptor = (uint8_t*)malloc(block_length * sizeof(data_type) * 4+1000);
        uint8_t* res = descriptor;
        res = ef->dump(descriptor);
        uint32_t segment_size = res - descriptor;
        descriptor = (uint8_t*)realloc(descriptor, segment_size);
        block_start_vec.push_back(descriptor);
        totalsize += segment_size;

        // elias_fanos.emplace_back(ef);
        // totalsize +=succinct::mapper::size_tree_of(*ef)->size;
    
    }
    
    double end_cr = getNow();
    double compress_time = end_cr - start_cr;
    double compress_throughput = N*sizeof(data_type) / (compress_time*1000000000);


    double compressrate = (totalsize)*100.0 / (sizeof(data_type) * N * 1.0);

    std::vector<succinct::elias_fano*>  elias_fanos;
    double decode_all_time = 0;
    double start = getNow();
    for(int i = 0; i< blocks;i++){
        succinct::elias_fano* ef = new succinct::elias_fano();
        ef->rebuild(block_start_vec[i]);
        elias_fanos.push_back(ef);
    }
    uint64_t mark_da = 0;
    block_length = block_size;
    for(int i=0;i<blocks;i++){
        // succinct::elias_fano ef = *(elias_fano_builders[i]);
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        succinct::elias_fano::select_enumerator it(*elias_fanos[i], 0);
        for (size_t i = 0; i < block_length; ++i) {
            mark_da+=it.next();
        }
    }
    
    double end = getNow();
    std::ofstream outfile("fix_log", std::ios::app);
    outfile<<mark_da<<std::endl;
    elias_fanos.clear();
    decode_all_time = end - start;
    double da_ns = decode_all_time / N * 1000000000;

    bool flag = true;
    std::vector<data_type> recover(data.size());
    double totaltime = 0.0; 

    // std::cout << "random access decompress!" << std::endl;
    std::vector<data_type> buffer(data.size());
    std::vector<uint32_t> ra_pos;
    for(int i=0;i<N;i++)
    {
        ra_pos.push_back(random(N));
        // ra_pos.push_back(i);
    }

    double randomaccesstime = 0.0;
    start = getNow();
    uint32_t mark = 0;
    int search_count = N;
    for(int i = 0; i< blocks;i++){
        succinct::elias_fano* ef = new succinct::elias_fano();
        ef->rebuild(block_start_vec[i]);
        elias_fanos.push_back(ef);
    }
    for (auto index: ra_pos)
    {

        // succinct::elias_fano ef;
        // ef.rebuild(block_start_vec[index/block_size]);
        // data_type tmpvalue = ef.select(index%block_size);
        data_type tmpvalue = elias_fanos[index/block_size]->select(index%block_size);

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
    double ra_ns = randomaccesstime / N * 1000000000;
    outfile<<mark<<std::endl;

    std::cout<<"Elias-Fano"<<" "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<0<<" "<<compressrate<<" "<<da_ns<<" "<<ra_ns<<" "<<compress_throughput<<std::endl;


}
