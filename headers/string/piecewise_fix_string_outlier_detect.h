
#ifndef PIECEWISEFIX_STRING_OUTLIER_H_
#define PIECEWISEFIX_STRING_OUTLIER_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string.h"
#include "string_utils.h"
#include "bit_read_string.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_int/serialize.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#define INF 0x7f7fffff

namespace Codecset
{
    int cal_bits(long_int delta, int max_string)
    {
        int tmp_bit = bits_long(delta) + 1;
        tmp_bit = std::min(max_string, tmp_bit);
        return tmp_bit;
    }
    // bit + theta0 + theta1 + delta
    // max_string is the max string size (each string in string_vec is padding to this size, calculate by bit)
    uint8_t *encodeArray8_string(std::vector<std::string> &string_vec, int start_idx, const size_t length, uint8_t *res, size_t max_string /*, long_int & theta0, long_int & theta1*/)
    {
        uint8_t *out = res;
        std::vector<long_int> ascii_vec;
        std::vector<long_int> index;

        for (int i = 0; i < length; i++)
        {
            ascii_vec.emplace_back(convertToASCII(string_vec[i + start_idx]));
            index.emplace_back(long_int(i));
        }

        string_lr mylr;
        mylr.caltheta(index, ascii_vec, length);

        long_float theta0 = mylr.theta0;
        long_int theta0_int = theta0.convert_to<long_int>();
        long_float theta1 = mylr.theta1;
        long_int theta1_int = theta1.convert_to<long_int>();

        std::vector<long_int> delta;
        long_int max_delta = 0;

        // a counter to plot the histogram
        int counter[max_string] = {0};
        // because strings are not fixed length, we take 256 as a back-up.
        for (auto i = 0; i < length; i++)
        {
            long_int tmp_val = ascii_vec[i] - (theta1_int * i + theta0_int);
            // std::cout<<"delta "<<i<<" : "<<tmp_val<<std::endl;
            delta.emplace_back(tmp_val);
            counter[cal_bits(abs(tmp_val), max_string)]++;
            if (abs(delta[i]) > max_delta)
            {
                max_delta = abs(delta[i]);
            }
        }

        int quantile_sum = 0;
        int threshold = 0;
        int compress_len = 0;
        int compress_min = length * max_string;
        int max_bit = 0;
        int temp = ceil((double)length / 64.);
        for (int i = 0; i < max_string; i++)
        {
            quantile_sum += counter[i];
            compress_len = quantile_sum * i + (length - quantile_sum) * max_string;
            if (quantile_sum == length)
            {
                compress_len = quantile_sum * i - temp * 12 * 8 - 12 * 8; // bitmap & lookup table
                max_bit = i;
            }
            if (compress_len < compress_min)
            {
                compress_min = compress_len;
                threshold = i;
            }
            if (quantile_sum == length)
            {
                break;
            }
        }

        // GONNA USE PIECEWISE, BECAUSE NOT A SINGLE OUTLIER OCCURS
        if (threshold == max_bit)
        {
            max_bit = bits_long(max_delta) + 1;
            memcpy(out, &max_bit, sizeof(uint32_t));
            out += sizeof(uint32_t);
            // models[nvalue*2] = theta0;
            // models[nvalue*2+1] = theta1;
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

            out = write_string_delta_string(delta.data(), out, max_bit, length);

            return out;
        }
        else
        {
            // GONNA USE OUTLIER DETECT
            uint64_t writebitmap[temp] = {0};
            bool vote[length];
            int usedData = 0;
            for (int i = 0; i < length; i++)
            {
                if (cal_bits(abs(delta[i]), max_string) > threshold)
                {
                    vote[i] = 0;
                }
                else
                {
                    vote[i] = 1;
                    usedData++;
                }
            }
            uint64_t tmpbit = 0;
            int wtitebittmp = 0;
            int outlier_num = length - usedData;
            bool have_outlier = true;
            std::vector<long_int> tmpoutlier;
            if (!outlier_num)
            {
                have_outlier = false;
            }

            max_delta = 0;
            std::vector<long_int> deltax;
            int wtiteoutlier = 0;
            for (int i = 0; i < length; i++)
            {

                if (vote[i])
                {
                    long_int tmp = delta[i];
                    deltax.emplace_back(tmp);
                    if (abs(tmp) > max_delta)
                    {
                        max_delta = abs(tmp);
                    }
                }
                else
                {
                    tmpbit += ((1L) << (63 - i % 64));
                    deltax.emplace_back((long_int)1);
                    tmpoutlier.emplace_back(convertToASCII(string_vec[i + start_idx]));
                }

                if (i % 64 == 63)
                {
                    writebitmap[wtitebittmp] = tmpbit;
                    wtitebittmp++;
                    tmpbit = 0;
                }
            }
            if (length % 64 != 0)
            {
                writebitmap[wtitebittmp] = tmpbit;
                wtitebittmp++;
            }
            int rank_lut_[temp] = {0};
            int basic_block_size_ = 64;
            initRankLut(length, writebitmap, rank_lut_, basic_block_size_);

            // WRITE IN OUTLIER FORM
            // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
            // ( 1L<<7 is to mark if we are using outlier or not)
            out[0] = (uint32_t)max_bit + (1L << 31);
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

            memcpy(out, &outlier_num, sizeof(outlier_num));
            out += sizeof(outlier_num);
            memcpy(out, &basic_block_size_, sizeof(basic_block_size_));
            out += sizeof(basic_block_size_);
            uint8_t *outlier_pos = out;
            out += sizeof(basic_block_size_);

            memcpy(out, writebitmap, temp * sizeof(writebitmap[0]));

            out += sizeof(writebitmap[0]) * temp;
            memcpy(out, rank_lut_, temp * sizeof(rank_lut_[0]));
            out += sizeof(rank_lut_[0]) * temp;

            out = write_string_delta_string(deltax.data(), out, max_bit, length);

            int shift = out - res;

            memcpy(outlier_pos, &shift, sizeof(shift));
            if (have_outlier)
            {
                for (auto outlier : tmpoutlier)
                {
                    mpz_t z;
                    mpz_init(z);
                    mpz_set(z, outlier.backend().data());
                    mpz_export(out, NULL, -1, 1, 0, 0, z);
                    out += sizeof(long_int);
                    mpz_clear(z);
                }
            }

            return out;
        }
    }

    // max_string is the length of each string, measured in bit
    void decodeArray8_string(uint8_t *in, const size_t length, long_int *out, size_t max_string, std::vector<std::string> &string_vec)
    {
        // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
        uint32_t whether_outlier;
        uint8_t *tmpin = in;
        whether_outlier = (tmpin[0] >> 31); // 1000 0000
        if (whether_outlier)
        {
            uint32_t maxerror;
            int *outlier_num = 0;
            maxerror = tmpin[0] - (1L << 31);
            tmpin += sizeof(maxerror);
            std::cout << maxerror << std::endl;

            tmpin += sizeof(uint32_t);
            uint32_t theta0_len;
            memcpy(&theta0_len, tmpin, sizeof(uint32_t));
            tmpin += sizeof(uint32_t);
            // theta0 = reinterpret_cast<long_int*>(tmpin);
            mpz_t tmp;
            mpz_init(tmp);
            mpz_import(tmp, theta0_len, -1, 1, 0, 0, tmpin);
            long_int theta0(tmp);
            // std::cout<<theta0<<std::endl;

            tmpin += theta0_len;

            uint32_t theta1_len;
            memcpy(&theta1_len, tmpin, sizeof(uint32_t));
            tmpin += sizeof(uint32_t);

            mpz_import(tmp, theta1_len, -1, 1, 0, 0, tmpin);
            long_int theta1(tmp);
            // std::cout<<theta1<<std::endl;
            mpz_clear(tmp);
            tmpin += theta1_len;

            outlier_num = reinterpret_cast<int *>(tmpin);
            tmpin += 8;
            int *outlier_position;
            outlier_position = reinterpret_cast<int *>(tmpin);
            tmpin += 4;
            uint8_t *bitmap_pos = tmpin;
            int temp = ceil((double)length / 64.);
            tmpin += temp * 12;

            uint8_t *outlier_pos = in + outlier_position[0];
            read_outlier_detect_string(tmpin, 0, 0, length, maxerror, theta1, theta0, out, outlier_pos, outlier_num[0], bitmap_pos, max_string);

            for (int i = 0; i < length; i++)
            {
                std::string tmp_string = convertToString(out[0]);
                string_vec.emplace_back(tmp_string);
                out++;
            }
        }
        else
        {
            uint32_t maxerror;
            maxerror = tmpin[0];
            tmpin += sizeof(maxerror);

            tmpin += sizeof(uint32_t);
            uint32_t theta0_len;
            memcpy(&theta0_len, tmpin, sizeof(uint32_t));
            tmpin += sizeof(uint32_t);
            // theta0 = reinterpret_cast<long_int*>(tmpin);
            mpz_t tmp;
            mpz_init(tmp);
            mpz_import(tmp, theta0_len, -1, 1, 0, 0, tmpin);
            long_int theta0(tmp);
            // std::cout<<theta0<<std::endl;

            tmpin += theta0_len;

            uint32_t theta1_len;
            memcpy(&theta1_len, tmpin, sizeof(uint32_t));
            tmpin += sizeof(uint32_t);

            mpz_import(tmp, theta1_len, -1, 1, 0, 0, tmpin);
            long_int theta1(tmp);
            // std::cout<<theta1<<std::endl;
            mpz_clear(tmp);
            tmpin += theta1_len;

            read_all_bit_fix_string(tmpin, 0, 0, length, maxerror, theta1, theta0, out);

            for (int i = 0; i < length; i++)
            {
                std::string tmp_string = convertToString(out[0]);
                string_vec.emplace_back(tmp_string);
                out++;
            }
        }
    }

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
