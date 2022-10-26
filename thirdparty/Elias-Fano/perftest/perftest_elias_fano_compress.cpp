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


int main(int argc, const char* argv[])
{
    std::string source_file = std::string(argv[1]);
    // We pick a CODEC
    std::vector<data_type> data = load_data_binary<data_type>("../integer_data/"+source_file);

    // std::vector<uint32_t> data;
    // std::ifstream srcFile("../integer_data/"+source_file, std::ios::in);
    // if (!srcFile.is_open())
    // {
    //     std::cout << "error opening source file." << std::endl;
    //     return 0;
    // }
    // while (srcFile.good())
    // {

    //     uint32_t next;
    //     srcFile >> next;
    //     // std::cout<<next<<std::endl;
    //     if (!srcFile.good())
    //     {
    //         break;
    //     }
    //     data.push_back(next);
    // }
    // srcFile.close();

    int N = data.size();
    if (data.size() == 0)
    {
        std::cout << "Empty vector" << std::endl;
        return 0;
    }
    // std::cout << "vector size = " << data.size() << std::endl;
    // std::cout << "vector size = " << data.size() * sizeof(uint32_t) / 1024.0 << "KB"
    //           << std::endl;

    int blocks = 1;
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size
    // std::cout << "Total blocks " << blocks << " block size " << block_size << std::endl;
    int block_length = block_size;
    uint64_t max_element = data[N-1];
    // std::vector<succinct::elias_fano *> elias_fano_builders;
    uint64_t totalsize = 0;
    // for(int i=0;i<blocks;i++)
    // {
    //     if (i == blocks - 1)
    //     {
    //         block_length = N - (blocks - 1) * block_size;
    //     }
    //     succinct::elias_fano::elias_fano_builder bvb(max_element, block_length);
    //     for(int j=0;j<block_length;j++)
    //     {
    //         bvb.push_back(data[i*block_size+j]);
    //     } 
    //     succinct::elias_fano ef(&bvb);
    //     totalsize +=succinct::mapper::size_tree_of(ef)->size;
    //     std::cout<<totalsize<<std::endl;
    //     elias_fano_builders.push_back(new succinct::elias_fano(&bvb));
    // }
    double start_cr = getNow();
    succinct::elias_fano::elias_fano_builder bvb(max_element, block_length);
    for(int j=0;j<N;j++)
    {
        bvb.push_back(data[j]);
    } 
    succinct::elias_fano ef(&bvb);
    totalsize +=succinct::mapper::size_tree_of(ef)->size;
    double end_cr = getNow();
    double compress_time = end_cr - start_cr;
    double compress_throughput = N*sizeof(data_type) / (compress_time*1000000000);


    double compressrate = (totalsize)*100.0 / (sizeof(data_type) * N * 1.0);
    // std::cout<<"compressed size "<<totalsize<<" uncompressed size "<<4*N<<std::endl;
    // std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;
    double decode_all_time = 0;
    double start = getNow();
    succinct::elias_fano::select_enumerator it(ef, 0);
    uint64_t mark_da = 0;
	for (size_t i = 0; i < N; ++i) {
	    mark_da+=it.next();
	}
    double end = getNow();

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
    }

    double randomaccesstime = 0.0;
    start = getNow();
    uint32_t mark = 0;
    int search_count = N;
    for (auto index: ra_pos)
    {
        // std::cout<<i<<std::endl;
        
        // int index = i;
        // uint32_t tmpvalue = (*elias_fano_builders[index/block_size]).select(index%block_size);
        uint32_t tmpvalue = ef.select(index);

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

    std::cout<<"Elias-Fano"<<" "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<0<<" "<<compressrate<<" "<<da_ns<<" "<<ra_ns<<" "<<compress_throughput<<std::endl;


    // std::cout << "random decoding time per int: " << std::setprecision(8)
    //           << randomaccesstime / search_count * 1000000000 << " ns" << std::endl;
    // std::cout << "random decoding speed: " << std::setprecision(10)
    //           << search_count / (randomaccesstime * 1000) << std::endl;

}
