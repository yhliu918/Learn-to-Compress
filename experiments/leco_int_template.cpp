#include "common.h"
#include "codecfactory.h"
#include "caltime.h"
#include "lr.h"
#include "piecewise_fix_integer_template.h"
#include "piecewise_fix_integer_template_float.h"
#include "piecewise_cost_integer_template.h"
#include "FOR_integer_template.h"
#include "delta_integer_template.h"

typedef uint64_t leco_type;

int random(int m)
{
  return rand() % m;
}

template <typename T>
static std::vector<T> load_data(const std::string& filename,
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
    // std::ofstream out("/home/lyh/Learn-to-Compress/integer_data/wiki_200M_uint64.txt",std::ios::out);
    // std::ios::sync_with_stdio(false);
    // for (auto i = 0; i < size; i++) {
    //   out << data[i] << std::endl;
    // }


  return data;
}

int main()
{
    using namespace Codecset;
    Leco_int<leco_type> codec;
    // alternatives : Leco_int, Delta_int, FOR_int

    std::vector<leco_type> data = load_data<leco_type>("/home/lyh/Learn-to-Compress/integer_data/wiki_200M_uint64", true);
    int N = data.size();

    std::cout << "vector size = " << data.size() << std::endl;
    std::cout << "vector size = " << data.size() * sizeof(leco_type) / 1024.0 << "KB"
              << std::endl;

    int blocks = 1000000;
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size

    // if using auto segmentation codecs
    // int delta = (1<<11);
    // codec.init(blocks, block_size, delta);

    std::cout << "Total blocks " << blocks << " block size " << block_size << std::endl;

    std::vector<uint8_t *> block_start_vec;

    uint64_t totalsize = 0;
    for (int i = 0; i < blocks; i++)
    {
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(leco_type)*4);
        uint8_t *res = descriptor;

        // std::cout<<data[i*block_size]<<" "<<data[i * block_size+block_size-1]<<std::endl;
        res = codec.encodeArray8_int(data.data() + (i * block_size), block_length, descriptor, i);
        uint32_t segment_size = res - descriptor;
        descriptor = (uint8_t *)realloc(descriptor, segment_size);
        block_start_vec.push_back(descriptor);
        totalsize += segment_size;
    }
    // if(totalsize == 0){
    //     totalsize = codec.get_block_nums();
    // }
    double compressrate = (totalsize)*100.0 / (sizeof(leco_type) * N * 1.0);
    std::cout<<"compressed size "<<totalsize<<" uncompressed size "<<sizeof(leco_type)*N<<std::endl;
    std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;
    bool flag = true;
    double totaltime = 0.0; 

    std::cout << "random access decompress!" << std::endl;
    double randomaccesstime = 0.0;
    double start = getNow();
    leco_type mark = 0;
    int search_count = N;
    for (int i = 0; i < search_count; i++)
    {

        // int index = random(N);
        int index = i;

        leco_type tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)index / block_size], index % block_size, NULL, N);

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
    double end = getNow();
    randomaccesstime += (end - start);

    std::cout << "random decoding time per int: " << std::setprecision(8)
              << randomaccesstime / search_count * 1000000000 << " ns" << std::endl;
    std::cout << "random decoding speed: " << std::setprecision(10)
              << search_count / (randomaccesstime * 1000) << std::endl;

    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}