
#include "codecfactory.h"
#include <cassert>
using namespace Codecset;

int randomint(int m)
{
    return rand() % m;
}

struct codec_vote{
    int select;
    double compression_rate;
};
codec_vote pick_block_size(int bsize [],int choice, int times, double sample_rate,int N, uint32_t * data, std::string codec_name){
    
    int select_= -1;

    int delta = 0;
    double mini_size = 1.0;
    std::vector<int> pick_vec;
    double sample_size_percent = sample_rate;
    int sample_size = N * sample_size_percent;

    for (int p = 0; p < times; p++)
    {
        pick_vec.push_back(randomint(1 / sample_size_percent));
    }
    for(int j = 0; j < choice; ++j)
    {
        IntegerCODEC &codec = *CODECFactory::getFromName(codec_name);
        int block_size = bsize[j];
        int blocks = N / block_size;
        codec.init(blocks, block_size, delta);
        int totalsize = 0;
        
        for (int p = 0; p < times; p++)
        {
            int pick_block = pick_vec[p];
            for (int i = 0; i < sample_size / block_size; i++)
            {
                uint8_t *descriptor = (uint8_t *)malloc(block_size * sizeof(uint64_t)+1024);
                uint8_t *res = descriptor;
                res = codec.encodeArray8(data + pick_block * sample_size + (i * block_size), block_size, descriptor, i);
                totalsize += (res - descriptor);
                free(descriptor);
            }
            if (sample_size % block_size != 0)
            {
                int block_length = sample_size - (sample_size / block_size) * block_size;
                uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t)+1024);
                uint8_t *res = descriptor;
                res = codec.encodeArray8(data+ (sample_size / block_size) * block_size, block_length, descriptor, sample_size / block_size);
                totalsize += (res - descriptor);
                free(descriptor);
            }
        }
        double compressrate = totalsize * 1.0 / (4 * times * sample_size * 1.0);
        if (compressrate < mini_size)
        {
            mini_size = compressrate;
            select_ = j;
        }

    }
    assert (select_ != -1);
    codec_vote tmp;
    tmp.select = select_;
    tmp.compression_rate = mini_size;
    return tmp;
}


