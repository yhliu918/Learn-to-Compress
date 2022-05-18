
#ifndef PIECEWISEFIX_STRING_TEMPLATE_H_
#define PIECEWISEFIX_STRING_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff
#include <limits>

namespace Codecset
{
    template <typename T>
    class Leco_integer
    {
    public:
    uint8_t* encodeArray8(T* in, const size_t length, uint8_t* res,
                        size_t nvalue) {
    if(length == 1){
      memcpy(res, in, sizeof(T));
      return res + sizeof(T);
    }
    int* indexes = new int[length];
    T* keys = new T[length];
    std::vector<bool> signvec;
    // double *keys_sample = new double [length];
    // double *indexes_sample = new double[length];
    uint8_t* out = res;
    for (uint32_t i = 0; i < length; i++) {
      indexes[i] = i;
      keys[i] = in[i];
    }
    T* delta = new T[length];

    lr_T<T> mylr;
    // mylr.caltheta_LOO(indexes,keys,length);
    mylr.caltheta(indexes, keys, length);

    delete[] indexes;
    delete[] keys;
    T max_error = 0;
    for (int i = 0; i < length; i++) {
        T pred = mylr.theta0 + mylr.theta1 * i;
        T tmp_val;
      if (pred > keys[i]) {
        tmp_val = pred - keys[i];
        signvec.emplace_back(false);
      } else if (pred < keys[i]) {
        tmp_val = keys[i] - pred;
        signvec.emplace_back(true);
      } else {
        tmp_val = 0;
        signvec.emplace_back(false);
      }
      delta[i] = tmp_val;

      if (tmp_val > max_error) {
        max_delta = tmp_val;
      }
    }
    // std::cout<<"delta[1] "<<delta[1]<<std::endl;
    uint32_t max_bit = bits_T<T>(max_error) + 1;
    memcpy(out, &max_bit, sizeof(uint32_t));
    out += sizeof(uint32_t);
    memcpy(out, &theta0, sizeof(T));
    out += sizeof(T);
    memcpy(out, &theta1, sizeof(T));
    out += sizeof(T);

    out = write_delta_T<T>(delta,signvec, out, max_bit, length);

    delete[]delta;

    return out;
  }



        uint8_t randomdecodeArray8_string(const char* in, int idx, T* result, int length) {
            const uint8_t* tmpin = reinterpret_cast<const uint8_t*>(in);
            uint32_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint32_t));
            tmpin += sizeof(uint32_t);
            // std::cout<<"max bit "<<maxbits<<std::endl;
            uint8_t block_pad_length;
            memcpy(&block_pad_length, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);

            T theta0;
            memcpy(&theta0, tmpin, sizeof(T));
            tmpin += sizeof(T);

            T theta1;
            memcpy(&theta1, tmpin, sizeof(T));
            tmpin += sizeof(T);

            uint8_t ori_length = 0;
            read_bit_fix_string<T>(tmpin, maxbits, idx, theta1, theta0, result, &ori_length);

            int shift = (block_pad_length - ori_length);
            *result = *result >> (uint8_t)(8 * shift);
            *result = *result << (uint8_t)(8 * shift);
            int tmp_length = block_pad_length - length;
            if (tmp_length > 0) {
                *result = *result >> (uint8_t)(8 * tmp_length);
            }
            else {
                *result = *result << (uint8_t)(8 * (-tmp_length));
            }
            //*result = *result >> (uint8_t)(8 * std::max(block_pad_length - length, 0));

            return ori_length;
        }

        /*
                void randomdecodeArray8_bsearch(int block_ind, int idx, T* result, int length)
                {

                    uint8_t* tmpin = block_start_vec[block_ind];
                    uint32_t maxbits;
                    memcpy(&maxbits, tmpin, sizeof(uint32_t));
                    tmpin += sizeof(uint32_t);
                    // std::cout<<"max bit "<<maxbits<<std::endl;

                    uint8_t block_pad_length;
                    memcpy(&block_pad_length, tmpin, sizeof(uint8_t));
                    tmpin += sizeof(uint8_t);

                    T theta0;
                    memcpy(&theta0, tmpin, sizeof(T));
                    tmpin += sizeof(T);

                    T theta1;
                    memcpy(&theta1, tmpin, sizeof(T));
                    tmpin += sizeof(T);

                    uint8_t ori_length = 0;
                    read_bit_fix_string<T>(tmpin, maxbits, idx, theta1, theta0, result, &ori_length);


                    int shift = (block_pad_length - ori_length);
                    *result = *result >> (8 * shift);
                    for (int i = 0; i < shift; i++) {
                        *result = (*result << 8) + (uint8_t)2;
                    }

                    *result = *result >> (8 * (block_pad_length - length));
                    // need to put the last few byte to padding_char



                    return;

                }
        */


        void setBlockSize(int block_size, int blocks)
        {
            this->block_size = block_size;
            this->blocks = blocks;
            return;
        }

        void push_back_block(uint8_t* block)
        {
            block_start_vec.emplace_back(block);
            return;
        }

        void push_back_firstkey(std::string& firstkey)
        {
            firstkey_each_block.append(firstkey);
            return;
        }

        void get_uncompressed_size(int& withpad, int& withoutpad)
        {
            withpad = totalsize_with_padding;
            withoutpad = totalsize_without_padding;
            std::cout << "totalsize_with_padding " << totalsize_with_padding << std::endl;
            std::cout << "totalsize_without_padding " << totalsize_without_padding << std::endl;

        }
        int get_max_padding_length()
        {
            return max_padding_length;
        }


    private:
        std::string firstkey_each_block;
        std::vector<uint8_t*> block_start_vec;
        int block_size;
        int blocks;
        int totalsize_without_padding;
        int totalsize_with_padding;
        std::vector<int> padding_length;
        std::vector<uint8_t> string_wo_padding_length;
        std::vector<std::string> padding_max;
        std::vector<std::string> padding_min;
        int max_padding_length;
        char* data_;
    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
