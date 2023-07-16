#include "common.h"
#include "codecfactory.h"
#include "caltime.h"
#include "lr.h"
#include "piecewise_fix_integer_template.h"
#include "piecewise_fix_integer_template_float.h"
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

int main(int argc, const char* argv[])
{
    using namespace Codecset;
    Leco_int<leco_type> codec;
    std::string method = "uncompressed";
    std::string source_file = std::string(argv[1]);
    int block_size = atoi(argv[2]);
    bool binary = atoi(argv[3]);
    std::string dictionary_file = std::string(argv[4]);
    // alternatives : Delta_int, Delta_cost, Delta_cost_merge, FOR_int, Leco_int, Leco_cost, Leco_cost_merge_hc,  Leco_cost_merge, Leco_cost_merge_double

    // load dictionary index 
    std::vector<leco_type> index;
    if(binary){
        index = load_data_binary<leco_type>(source_file);
    }
    else{
        index = load_data<leco_type>(source_file);
    }

    // load dictionary
    std::vector<leco_type> dictionary;
    dictionary = load_data<leco_type>(dictionary_file);
    int dict_size = dictionary.size();
    int N = index.size();
    int blocks = dict_size / block_size;
    if (blocks * block_size < dict_size)
    {
        blocks++;
    } // handle with the last block, maybe < block_size
    codec.init(blocks, block_size);

    std::vector<uint8_t*> block_start_vec;

    


    double origin_size = (sizeof(leco_type) * dict_size * 1.0);
    double compressrate = 1;

    double start = getNow();
    leco_type mark = 0;
    // lookup dictionary
    bool flag = true;
    int repeat = 1;
    for(auto tmp_index : index){
        leco_type tmpvalue = dictionary[tmp_index];
        mark += tmpvalue;

    }


    double end = getNow();
    double randomaccesstime = (end - start);
    std::ofstream outfile("fix_log", std::ios::app);
    outfile<<mark<<std::endl;

    double ra_ns = randomaccesstime / (N * repeat) * 1000000000;

    std::cout << method << " " << source_file << " " << blocks << " " << compressrate <<" " << ra_ns<< std::endl;


    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
