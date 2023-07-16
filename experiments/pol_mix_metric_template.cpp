#include "common.h"
#include "codecfactory.h"
#include "caltime.h"
#include "lr.h"
#include "poly_fix_integer_template.h"
#include "piecewise_fix_integer_template.h"
#include "piecewise_cost_merge_integer_template_double.h"
#include "FOR_integer_template.h"
#include "delta_integer_template.h"
#include "delta_cost_integer_template.h"
#include "delta_cost_merge_integer_template.h"
#include "piecewise_cost_merge_integer_template_test.h"

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
double linearity(std::vector<int64_t> data, int64_t max_val, int64_t min_val){
    lr_int_T<int64_t> mylr;
    mylr.caltheta(data.data(), data.size());
    int64_t data_range = max_val - min_val;
    if (data_range == 0){
        return 0;
    }
    double metric = 0;
    for(int i = 0; i < data.size(); i++){
        double tmp_val =data[i] - (mylr.theta0 + mylr.theta1 * (double)i);
        metric += (double)abs(tmp_val) / (double)data_range;
    }
    return 1 - metric / (double)data.size();

}

double var(std::vector<int64_t> data, int64_t max_val, int64_t min_val, int64_t sum_val){
    double mid = (double)sum_val / (double)data.size();
    double metric = 0;
    int64_t data_range = max_val - min_val;
    if (data_range == 0){
        return 1;
    }
    for(int i = 0; i < data.size(); i++){
        double tmp_val = ((double)data[i] - mid);
        // std::cout<<data[i]<<" "<<mid<<std::endl;
        metric += (double)abs(tmp_val) / (double)data_range;
    }
    return 1 - metric / (double)data.size();

}

int main(int argc, const char* argv[])
{
    using namespace Codecset;
    
    std::string method = "Poly_mix";
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int model_size = atoi(argv[3]);
    bool binary = atoi(argv[4]);
    double threshold = atof(argv[5]);
    
    std::vector<leco_type> data;
    if(binary){
        data = load_data_binary<leco_type>(source_file);
        // data = load_data_binary<leco_type>("../data/" + source_file);
    }
    else{
        data = load_data<leco_type>(source_file);
        // data = load_data<leco_type>("../data/" + source_file);
    }
    int N = data.size();

    int block_size = data.size() / blocks;
    blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } // handle with the last block, maybe < block_size
    // convert block_size to const
    // const int contsblock_size = const_cast<const int&>(block_size);

    Leco_int_poly<leco_type, 2> codec_degree2;
    codec_degree2.init(blocks, block_size);

    Leco_int<leco_type> codec;
    codec.init(blocks, block_size);

    std::vector<uint8_t*> block_start_vec;
    int poly_degree_count[3] = {0, 0, 0};
    uint64_t totalsize = 0;
    for (int i = 0; i < blocks; i++)
    {
        // std::cout<<"block "<<i<<std::endl;
        int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
        std::vector<int64_t> delta_seq;
        int64_t max_val = INT64_MIN;
        int64_t min_val = INT64_MAX;
        int64_t sum_val = 0;
        for(int j=0; j< block_length-1;j++){
            int64_t tmp_delta = data[i * block_size + j+1] - data[i * block_size + j];
            delta_seq.push_back(tmp_delta);
            if(tmp_delta > max_val){
                max_val = tmp_delta;
            }
            if(tmp_delta < min_val){
                min_val = tmp_delta;
            }
            sum_val += tmp_delta;
        }
        int poly_degree = 1;
        double linearity_metric = linearity(delta_seq, max_val, min_val);
        double var_metric = var(delta_seq, max_val, min_val, sum_val);
        // std::cout<<linearity_metric<<" "<<var_metric<<std::endl;
        if (linearity_metric - var_metric> threshold){
            poly_degree = 2;
        }

        uint8_t* descriptor = (uint8_t*)malloc(block_length * sizeof(leco_type) * 4);
        uint8_t* res = descriptor;

    
        uint32_t segment_size = 0;
        if(poly_degree == 1){
            res = codec.encodeArray8_int(data.data() + (i * block_size), block_length, descriptor, i);
        }
        else{
            res = codec_degree2.encodeArray8_int(data.data() + (i * block_size), block_length, descriptor, i);
        }
        // std::cout<<poly_degree-1<<std::endl;
        poly_degree_count[poly_degree - 1] += 1;
        segment_size = res - descriptor;
        descriptor = (uint8_t*)realloc(descriptor, segment_size);
        block_start_vec.push_back(descriptor);
        totalsize += segment_size;
    }


    double origin_size = (sizeof(leco_type) * N * 1.0);
    double total_model_size = model_size * blocks;
    double cr_wo_model = (totalsize - total_model_size) * 100.0 / origin_size;
    double cr_model = total_model_size * 100.0 / origin_size;
    double compressrate = (totalsize) * 100.0 / origin_size;
    std::cout<< (double)poly_degree_count[0]/blocks<<" "<<(double)poly_degree_count[1]/blocks<<" "<<(double)poly_degree_count[2]/blocks<<std::endl;
    std::cout<<compressrate<<std::endl;
    // bool flag = true;
    // std::vector<leco_type> recover(data.size());
    // double totaltime = 0.0;
    // // std::cout << "decompress all!" << std::endl;
    // int repeat = 10;

    // double start = getNow();
    // for (int iter = 0; iter < repeat; iter++)
    // {
    //     for (int i = 0; i < blocks; i++)
    //     {
    //         int block_length = block_size;
    //         if (i == blocks - 1)
    //         {
    //             block_length = N - (blocks - 1) * block_size;
    //         }
    //         codec.decodeArray8(block_start_vec[i], block_length, recover.data() + i * block_size, i);
    //     }
    //     for (auto index : codec.mul_add_diff_set)
    //     {
    //         recover[index.first] += index.second;
    //     }
    //     #ifndef NDEBUG
    //     for (int j = 0; j < N; j++)
    //     {
    //         if (data[j ] != recover[j ])
    //         {
    //             std::cout<< " num: " << j << " true is: " << data[j] << " predict is: " << recover[j] << std::endl;
    //             std::cout << "something wrong! decompress all failed" << std::endl;
    //             flag = false;
    //             break;
    //         }
    //     }
    //     if (!flag)
    //     {
    //         break;
    //     }
    //     #endif
    // }
    // double end = getNow();
    // totaltime += (end - start);
    // double da_ns = totaltime / (N*repeat) * 1000000000;

    // // std::cout << "random access decompress!" << std::endl;
    // std::vector<uint32_t> ra_pos;
    // repeat = 1;
    // for (int i = 0;i < N * repeat;i++) {
    //     ra_pos.push_back(random(N));
    // }
    // flag = true;
    // double randomaccesstime = 0.0;
    // start = getNow();
    // leco_type mark = 0;
    // for (auto index : ra_pos)
    // {

    //     leco_type tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)index / block_size], index % block_size, NULL, N);
    //     mark += tmpvalue;
    //     #ifndef NDEBUG
    //     if (data[index] != tmpvalue)
    //     {

    //         std::cout << "num: " << index << "true is: " << data[index] << " predict is: " << tmpvalue << std::endl;
    //         flag = false;
    //         std::cout << "something wrong! random access failed" << std::endl;
    //     }
    //     if (!flag)
    //     {
    //         break;
    //     }
    //     #endif
    // }
    // end = getNow();
    // randomaccesstime += (end - start);
    // std::ofstream outfile("fix_log", std::ios::app);
    // outfile<<mark<<std::endl;

    // double ra_ns = randomaccesstime / (N * repeat) * 1000000000;

    // std::cout << method << " " << source_file << " " << blocks << " " << compressrate << " " << cr_model << " " << cr_wo_model << " " << da_ns << " " << ra_ns<<" "<<0<< std::endl;


    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
