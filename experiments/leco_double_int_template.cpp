#include "common.h"
#include "codecfactory.h"
#include "caltime.h"
#include "lr.h"
#include "piecewise_fix_integer_template.h"
#include "piecewise_fix_integer_template_float.h"
#include "piecewise_cost_integer_template.h"
#include "piecewise_cost_merge_integer_template_double.h"
#include "piecewise_cost_merge_integer_template_double_link.h"
#include "FOR_integer_template.h"
#include "delta_integer_template.h"
#include "delta_cost_integer_template.h"
#include "delta_cost_merge_integer_template.h"

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

int main(int argc, const char* argv[])
{
    using namespace Codecset;
    Leco_cost_merge_double_link<leco_type> codec;
    std::string method = "leco_cost";
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int delta = atoi(argv[3]);
    int model_size = atoi(argv[4]);
    // alternatives : Delta_int, Delta_cost, Delta_cost_merge, FOR_int, Leco_int, Leco_cost, Leco_cost_merge_hc,  Leco_cost_merge, Leco_cost_merge_double

    std::vector<leco_type> data = load_data<leco_type>("../data/" + source_file);

    int N = data.size();
    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size

    // if using auto segmentation codecs
    // int delta = 32;
    codec.init(blocks, block_size, delta);


    std::vector<uint8_t*> block_start_vec;
    double start_cr = getNow();
    uint64_t totalsize = 0;
    for (int i = 0; i < blocks; i++)
    {
        // std::cout<<"block "<<i<<std::endl;
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        uint8_t* descriptor = (uint8_t*)malloc(block_length * sizeof(leco_type) * 4);
        uint8_t* res = descriptor;
        // if adaptive segment
        res = codec.encodeArray8_int(data.data(), block_length, descriptor, i);
        // if fixed length segment
        // res = codec.encodeArray8_int(data.data()+(i*block_size), block_length, descriptor, i);
        uint32_t segment_size = res - descriptor;
        descriptor = (uint8_t*)realloc(descriptor, segment_size);
        block_start_vec.push_back(descriptor);
        totalsize += segment_size;
    }
    if (totalsize == 0) {
        totalsize = codec.get_total_byte();
    }
    double end_cr = getNow();
    double cr_through = N * sizeof(leco_type) / ((end_cr-start_cr) * 1000000000);


    blocks = codec.get_total_blocks();
    double origin_size = (sizeof(leco_type) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize) * 100.0 / origin_size;

    bool flag = true;
    std::vector<leco_type> recover(data.size());
    double totaltime = 0.0;
    // std::cout << "decompress all!" << std::endl;
    double start = getNow();
    codec.decodeArray8(N, recover.data(), N);
    // for (int j = 0; j < N; j++)
    // {
    //     if (data[j] != recover[j])
    //     {
    //         std::cout <<"num: " << j << " true is: " << data[j] << " predict is: " << recover[j] << std::endl;
    //         std::cout << "something wrong! decompress failed" << std::endl;
    //         flag = false;
    //         break;
    //     }
    // }
    double end = getNow();
    totaltime += (end - start);
    double da_ns = totaltime / N * 1000000000;

    // std::cout << "random access decompress!" << std::endl;
    int repeat = 1;
    std::vector<uint32_t> ra_pos;
    ra_pos.reserve(N*repeat);
    for(int i=0;i<N*repeat;i++){
        ra_pos.push_back(random(N));
        // ra_pos.push_back(i);
    }
    flag = true;
    double randomaccesstime = 0.0;
    codec.search_node.reserve(8);
    start = getNow();

    codec.art.Build(codec.art_build_vec);
    leco_type mark = 0;
    int segment_id = codec.get_segment_id(ra_pos[0]), next_segment_id = 0;

    for (int i=0;i<ra_pos.size();i++)
    {
        auto index=ra_pos[i];
        if(i<ra_pos.size()-1) next_segment_id = codec.get_segment_id(ra_pos[i+1]);

        // leco_type tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)index / block_size], index % block_size, NULL, N);
        leco_type tmpvalue = codec.randomdecodeArray8(segment_id, block_start_vec[(int)index / block_size], index, NULL, N);
        mark += tmpvalue;
        segment_id=next_segment_id;
        // if (data[index] != tmpvalue)
        // {
        //     std::cout << "num: " << index << "true is: " << data[index] << " predict is: " << tmpvalue << std::endl;
        //     std::cout << "something wrong! random access failed" << std::endl;
        //     flag = false;
        //     break;
        // }
    }
    
    end = getNow();
    randomaccesstime += (end - start);
    std::cout<< mark << std::endl;
    double ra_ns = randomaccesstime / (N*repeat) * 1000000000;

    
    std::cout<<method<<" "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<cr_model<<" "<<cr_wo_model<<" "<<da_ns<<" "<<ra_ns<<" "<<cr_through<<std::endl;


    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
