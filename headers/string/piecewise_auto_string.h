
#ifndef PIECEWISE_STRING_H_
#define PIECEWISE_STRING_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string.h"
#include "string_utils.h"
#include "bit_read_string.h"

namespace Codecset
{
    class Piecewise_auto
    {
    public:
        // start_index + bit + theta0 + theta1 + numbers + delta
        void init(long_int delta)
        {
            maxerror = delta;
        }

        uint32_t lower_bound(long_int v, uint32_t len)
        {
            uint32_t m;
            uint32_t x = 0;
            uint32_t y = len - 1;
            while (x <= y)
            {

                m = x + (y - x) / 2;
                if (v < segment_index[m])
                    y = m - 1;
                else
                    x = m + 1;
            }
            return y;
        }

        uint8_t *encodeArray8(std::vector<std::string> &string_vec, int start_idx, const size_t length, uint8_t *res, size_t nvalue)
        {

            std::vector<long_int> ascii_vec;
            std::vector<int> index;
            for (int i = 0; i < length; i++)
            {
                ascii_vec.emplace_back(convertToLongInt(string_vec[i + start_idx]));
                index.emplace_back(i);
            }

            long_int high_slope = inf;
            long_int low_slope = 0;
            long_int origin_key = ascii_vec[0];
            int origin_index = index[0];
            int end_index = index[0];
            int total_index = 0;
            for (int i = 1; i < length; i++)
            {
                long_int key = ascii_vec[i];
                int id = index[i];
                long_int tmp_point_slope = (key - origin_key) / (id - origin_index);
                //std::cout<<low_slope<<" "<<tmp_point_slope<<" "<<high_slope<<std::endl;

                if (tmp_point_slope >= low_slope && tmp_point_slope <= high_slope && (id - origin_index) <= 1000000)
                {
                    long_int tmp_high_slope = ((key + maxerror - origin_key)) / ((id - origin_index));
                    long_int tmp_low_slope = ((key - maxerror - origin_key)) / ((id - origin_index));

                    if (tmp_low_slope < 0)
                    {
                        tmp_low_slope = 0;
                    }
                    if (tmp_high_slope <= high_slope)
                    {
                        high_slope = tmp_high_slope;
                    }
                    if (low_slope <= tmp_low_slope)
                    {
                        low_slope = tmp_low_slope;
                    }
                    end_index = id;
                }
                else
                {

                    long_int slope = (high_slope + low_slope) / 2;
                    int max_error = 0;

                    if (end_index == origin_index)
                    {
                        slope = 1;
                    }
                    int seg_len = end_index - origin_index + 1;
                    long_int theta0_int = ascii_vec[origin_index];
                    long_int theta1_int = slope;

                    std::vector<long_int> delta;
                    long_int max_delta = 0;
                    for (int j = origin_index; j <= end_index; j++)
                    {
                        long_int tmp_val = ascii_vec[j] - (theta1_int * (j - origin_index) + theta0_int);
                        delta.emplace_back(tmp_val);
                        if (abs(tmp_val) > max_delta)
                        {
                            max_delta = abs(tmp_val);
                        }
                    }

                    uint32_t max_bit = 0;
                    if (max_delta)
                    {
                        max_bit = bits_long(max_delta) + 1;
                    }

                    uint8_t *descriptor = (uint8_t *)malloc((end_index - origin_index + 1) * sizeof(uint64_t) * 100);

                    uint8_t *out = descriptor;

                    memcpy(out, &origin_index, sizeof(int));
                    out += sizeof(int);

                    memcpy(out, &max_bit, sizeof(uint32_t));
                    out += sizeof(uint32_t);


                    mpz_t z;
                    mpz_init(z);
                    mpz_set(z, theta0_int.backend().data());
                    auto theta0_len = (mpz_sizeinbase(z, 2) + 7) / 8;
                    memcpy(out, &theta0_len, sizeof(uint32_t));
                    out += sizeof(uint32_t);
                    mpz_export(out, &theta0_len, -1, 1, 0, 0, z);
                    out += theta0_len;

                    mpz_set(z, theta1_int.backend().data());
                    auto theta1_len = (mpz_sizeinbase(z, 2) + 7) / 8;
                    memcpy(out, &theta1_len, sizeof(uint32_t));
                    out += sizeof(uint32_t);
                    mpz_export(out, &theta1_len, -1, 1, 0, 0, z);
                    out += theta1_len;
                    mpz_clear(z);

                    memcpy(out, &seg_len, sizeof(int));
                    out += sizeof(int);

                    out = write_string_delta_string(delta.data(), out, max_bit, seg_len);

                    descriptor = (uint8_t *)realloc(descriptor, (out - descriptor));
                    block_start_vec.push_back(descriptor);
                    segment_index.push_back(origin_index);

                    total_byte += (out - descriptor);

                    high_slope = inf;
                    low_slope = 0;
                    origin_index = id;
                    origin_key = key;
                    end_index = id;
                }
            }

            long_int slope = (high_slope + low_slope) / 2;
            int max_error = 0;

            if (end_index == origin_index)
            {
                slope = 1;
            }
            int seg_len = end_index - origin_index + 1;

            long_int theta0_int = ascii_vec[origin_index];
            long_int theta1_int = slope;

            std::vector<long_int> delta;
            long_int max_delta = 0;

            for (auto j = origin_index; j <= end_index; j++)
            {
                long_int tmp_val = ascii_vec[j] - (theta1_int * (j - origin_index) + theta0_int);
                delta.emplace_back(tmp_val);
                if (abs(tmp_val) > max_delta)
                {
                    max_delta = abs(tmp_val);
                }
            }

            uint32_t max_bit = 0;
            if (max_delta)
            {
                max_bit = bits_long(max_delta) + 1;
            }

            uint8_t *descriptor = (uint8_t *)malloc(seg_len * sizeof(uint64_t) * 100);
            uint8_t *out = descriptor;

            memcpy(out, &origin_index, sizeof(int));
            out += sizeof(int);

            memcpy(out, &max_bit, sizeof(uint32_t));
            out += sizeof(uint32_t);

            mpz_t z;
            mpz_init(z);
            mpz_set(z, theta0_int.backend().data());
            auto theta0_len = (mpz_sizeinbase(z, 2) + 7) / 8;
            memcpy(out, &theta0_len, sizeof(uint32_t));
            out += sizeof(uint32_t);
            mpz_export(out, &theta0_len, -1, 1, 0, 0, z);
            out += theta0_len;

            mpz_set(z, theta1_int.backend().data());
            auto theta1_len = (mpz_sizeinbase(z, 2) + 7) / 8;
            memcpy(out, &theta1_len, sizeof(uint32_t));
            out += sizeof(uint32_t);
            mpz_export(out, &theta1_len, -1, 1, 0, 0, z);
            out += theta1_len;
            mpz_clear(z);

            memcpy(out, &seg_len, sizeof(int));
            out += sizeof(int);

            out = write_string_delta_string(delta.data(), out, max_bit, seg_len);

            descriptor = (uint8_t *)realloc(descriptor, (out - descriptor));
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);

            total_byte += (out - descriptor);
            return res;
        }

        void decodeArray8(uint8_t *in, int length, long_int *out, size_t nvalue, std::vector<std::string> &string_vec)
        {
            // start_index + bit + theta0 + theta1 + numbers + delta
            long_int *tmpout = out;
            int len = block_start_vec.size();
            for (int i = 0; i < len; i++)
            {
                //if(i==len){return;}
                uint8_t *this_block = block_start_vec[i];
                uint8_t *tmpin = this_block;
                uint32_t maxbits;
                int start_ind;
                int numbers;
                
                memcpy(&start_ind, tmpin, sizeof(uint32_t));
                tmpin += sizeof(uint32_t);
                memcpy(&maxbits, tmpin, sizeof(uint32_t));
                tmpin += sizeof(uint32_t);

                uint32_t theta0_len;
                memcpy(&theta0_len, tmpin, sizeof(uint32_t));
                tmpin += sizeof(uint32_t);
                mpz_t tmp;
                mpz_init(tmp);
                mpz_import(tmp, theta0_len, -1, 1, 0, 0, tmpin);
                long_int theta0(tmp);
                //std::cout<<theta0<<std::endl;

                tmpin += theta0_len;

                uint32_t theta1_len;
                memcpy(&theta1_len, tmpin, sizeof(uint32_t));
                tmpin += sizeof(uint32_t);

                mpz_import(tmp, theta1_len, -1, 1, 0, 0, tmpin);
                long_int theta1(tmp);
                //std::cout<<theta1<<std::endl;
                mpz_clear(tmp);
                tmpin += theta1_len;

                memcpy(&numbers, tmpin, sizeof(int));
                tmpin += sizeof(int);
                //std::cout<<i<< " seg length "<<numbers<<std::endl;
                //std::cout<<i<<": [ "<<start_ind<<", "<<start_ind+numbers-1<<" ]"<<std::endl;

                if (numbers == 1)
                {
                    tmpout[0] = theta0;
                    tmpout++;
                }
                else
                {
                    //read_all_bit_fix_string(tmpin, 0, 0, numbers, maxbits, theta1, theta0, tmpout);
                    tmpout+=numbers;
                }
            }
            // for (int i = 0; i < length; i++)
            // {
            //     std::string tmp_string = convertToString(out[0]);
            //     string_vec.emplace_back(tmp_string);
            //     out++;
            // }

            
        }
        uint32_t get_total_byte(){
            return total_byte; 
        }
        int get_total_seg(){
            return segment_index.size();
        }

    private:
        long_int inf = ((long_int)1 << 100);
        std::vector<uint8_t *> block_start_vec;
        std::vector<int> segment_index;
        uint32_t total_byte = 0;
        long_int maxerror = (1U << 10) - 1;

    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
