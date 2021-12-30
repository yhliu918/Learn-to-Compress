#ifndef SEARCH_BLOCKSIZE_H_
#define SEARCH_BLOCKSIZE_H_

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"

namespace Codecset
{
    class search_bsize
    {
    public:
        int search_block_size(std::string codec_name, std::vector<uint32_t> &data, 
                                double sample_size_percent, int times)
        {
            int bsize[7] = {200, 400, 800, 1600, 3200, 6400, 10000};
            int select_;

            int delta = 0;
            int N = data.size();

            double mini_size = 1.0;
            double start = getNow();
            std::vector<int> pick_vec;
            int sample_size = N * sample_size_percent;

            for (int j = 0; j < 7; ++j)
            {
                IntegerCODEC &codec = *CODECFactory::getFromName(codec_name);
                int block_size = bsize[j];
                int blocks = N / block_size;
                if(block_size*blocks<N) ++blocks;
                codec.init(blocks, block_size, delta);

                int totalsize = 0;
                for (int p = 0; p < times; p++)
                {
                    int sample_blocks = sample_size / block_size;
                    if (sample_size % block_size != 0) ++sample_blocks;
                    std::set<int> s;
                    while (s.size() < sample_blocks) {
                        s.insert(rand() % (N / block_size));
                    }
                    auto it = s.begin();
                    for (int i = 0; i < sample_blocks; i++)
                    {
                        int block_length = block_size;
                        if(i == sample_blocks-1) block_length = sample_size - i * block_size;
                        uint8_t *descriptor = (uint8_t *)malloc(2 * block_length * sizeof(uint64_t));
                        uint8_t *res = descriptor;
                        int t = *it;
                        ++it;
                        res = codec.encodeArray8(data.data() + t * block_size, block_length, descriptor, i);
                        descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
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
            //std::cout << "ok" << std::endl;
            mini_size = 1.0;
            int select_delta = 0;
            for (int j = -100; j <= 100; j += 50)
            {
                IntegerCODEC &codec = *CODECFactory::getFromName(codec_name);
                int block_size = bsize[select_] + j;
                int blocks = N / block_size;
                if(block_size*blocks<N) ++blocks;
                codec.init(blocks, block_size, delta);
                int totalsize = 0;
                for (int p = 0; p < times; p++)
                {
                    int sample_blocks = sample_size / block_size;
                    if (sample_size % block_size != 0) ++sample_blocks;
                    std::set<int> s;
                    while (s.size() < sample_blocks) {
                        s.insert(rand() % (N / block_size));
                    }
                    auto it = s.begin();
                    for (int i = 0; i < sample_blocks; i++)
                    {
                        int block_length = block_size;
                        if(i == sample_blocks-1) block_length = sample_size - i * block_size;
                        uint8_t *descriptor = (uint8_t *)malloc(2 * block_length * sizeof(uint64_t));
                        //std::cout << "malloc" << std::endl;
                        uint8_t *res = descriptor;
                        int t = *it;
                        ++it;
                        res = codec.encodeArray8(data.data() + t * block_size, block_length, descriptor, i);
                        descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
                        totalsize += (res - descriptor);
                        free(descriptor);
                    }
                }
                double compressrate = totalsize * 1.0 / (4 * times * sample_size * 1.0);
                if (compressrate < mini_size)
                {
                    mini_size = compressrate;
                    select_delta = j;
                }
            }

            double end = getNow();
            double search_time = end - start;

            return bsize[select_] + select_delta;
        }
    };
}

#endif