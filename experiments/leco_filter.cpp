#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
int random(int m)
{
    return rand() % m;
}
int filter_range(uint8_t* in, const size_t length, int filter, uint32_t* out) { 
    // only consider > filter, return [return_value, end] for demo
    int counter = 0;
    double theta0;
    float theta1;
    uint8_t maxerror;
    uint8_t* tmpin = in;
    maxerror = tmpin[0];
    tmpin++;
    memcpy(&theta0, tmpin, sizeof(theta0));
    tmpin += sizeof(theta0);
    memcpy(&theta1, tmpin, sizeof(theta1));
    tmpin += sizeof(theta1);
    int delta_interval = 0;
    if(maxerror){
        delta_interval = (1<<(maxerror-1));
    }
    int thre = std::max((double)(filter +1 - delta_interval - theta0)/theta1,0.);
    if(thre >= length){
        return counter;
    }
    else{
        if (maxerror) {
            counter = read_all_bit_fix_range<uint32_t>(tmpin ,0, thre, length, maxerror,theta1,theta0, out, filter);
        }
        else{
            for (int i = thre;i < length;i++) {
                uint32_t pred =  (long long)(theta0 + theta1*i);
                if(pred > filter){
                    out[0]= pred;
                    out++;
                    counter++;
                }
            }
        }
    }
    return counter;
}

int main(int argc, const char *argv[])
{
    using namespace Codecset;
    std::string method = "piecewise_fix_op_max_test";
    std::string source_file = std::string(argv[1]);
    int blocks = atoi(argv[2]);
    int model_size = atoi(argv[3]);
    int filter = atoi(argv[4]);
    std::ios::sync_with_stdio(false);
    // We pick a CODEC
    IntegerCODEC &codec = *CODECFactory::getFromName(method);

    std::vector<uint32_t> data;
    std::ifstream srcFile("" + source_file, std::ios::in);
    if (!srcFile)
    {
        std::cout << "error opening source file." << std::endl;
        return 0;
    }
    while (srcFile.good())
    {
        uint32_t next;
        srcFile >> next;
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
        uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t));
        uint8_t *res = descriptor;

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
    int repeat = 1;
    int total_counter = 0;

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
            total_counter+=filter_range(block_start_vec[i], block_length, filter, recover.data()+total_counter);
        }
    }
    end = getNow();
    // std::cout<<total_counter<<std::endl;
    totaltime += (end - start);


    std::cout<<method<<" "<<source_file<<" "<<filter<<" "<<blocks<<" "<<compressrate<<" "<<totaltime<<" "<<total_counter<<std::endl;
    
    for (int i = 0; i < (int)block_start_vec.size(); i++)
    {
        free(block_start_vec[i]);
    }
}
