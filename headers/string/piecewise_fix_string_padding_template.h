
#ifndef PIECEWISEFIX_STRING_TEMPLATE_H_
#define PIECEWISEFIX_STRING_TEMPLATE_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string.h"
#include "string_utils.h"
#include "bit_read_string.h"
#define INF 0x7f7fffff
#include <limits>

namespace Codecset
{
    template <typename T>
    class Leco_string
    {
    public:
        void Padd_one_string(std::string& s, char x) {
            s.append(max_padding_length - s.size(), x);
        }

        void Padding_string(std::vector<std::string>& string_vec, int N, char x) {
            totalsize_without_padding = 0;
            totalsize_with_padding = 0;
            max_padding_length = 0;
            for (int i = 0; i < blocks; i++) {
                int block_length = block_size;
                if (i == blocks - 1) {
                    block_length = N - (blocks - 1) * block_size;
                }
                uint8_t max_length = 0;

                for (int j = 0; j < block_length; j++) {
                    totalsize_without_padding += string_vec[i * block_size + j].size();
                    if (string_vec[i * block_size + j].size() > max_length) {
                        max_length = string_vec[i * block_size + j].size();
                    }
                }
                if (max_padding_length < max_length) {
                    max_padding_length = max_length;
                }
                padding_length.emplace_back(max_length);
                totalsize_with_padding += max_length * block_length;
                for (int j = 0; j < block_length; j++) {
                    int temp_len = string_vec[i * block_size + j].size();
                    string_wo_padding_length.emplace_back((uint8_t)temp_len);
                    std::string tmp_str_max = string_vec[i * block_size + j];
                    std::string tmp_str_min = tmp_str_max;
                    padding_max.emplace_back(tmp_str_max.append(
                        max_length - temp_len, std::numeric_limits<char>::max()));
                    padding_min.emplace_back(tmp_str_min.append(max_length - temp_len, 1));
                    string_vec[i * block_size + j].append(max_length - temp_len, x);
                }
            }
        }

        uint8_t* encodeArray8_string(std::vector<std::string>& string_vec,
            int start_idx, const size_t length, uint8_t* res,
            size_t nvalue) {
            uint8_t* out = res;
            std::vector<T> ascii_vec;
            std::vector<T> ascii_vec_min;
            std::vector<T> ascii_vec_max;
            std::vector<long_int>
                long_int_vec;  // because the calculation in linear regression may
                               // overflow, we choose to adopt long_int in the string_lr
                               // and convert its param to type T later
            std::vector<int> index;
            for (size_t i = 0; i < length; i++) {
                ascii_vec.emplace_back(convertToASCII<T>(string_vec[i + start_idx]));
                ascii_vec_min.emplace_back(convertToASCII<T>(padding_min[i + start_idx]));
                ascii_vec_max.emplace_back(convertToASCII<T>(padding_max[i + start_idx]));
                long_int_vec.emplace_back(convertToLongInt(string_vec[i + start_idx]));
                index.emplace_back(i);
            }

            string_lr mylr;
            mylr.caltheta(index, long_int_vec, length);

            T theta0 = 0;
            // auto theta0_len = (mpz_sizeinbase(mylr.theta0.backend().data(), 2) + 7) /
            // 8;
            mpz_export(&theta0, NULL, -1, 1, 1, 0, mylr.theta0.backend().data());
            // print_u128_u(theta0);
            // printf("\n");

            T theta1 = 0;
            // auto theta1_len = (mpz_sizeinbase(mylr.theta1.backend().data(), 2) + 7) /
            // 8;
            mpz_export(&theta1, NULL, -1, 1, 1, 0, mylr.theta1.backend().data());
            // print_u128_u(theta1);
            // printf("\n");
            // std::cout<<"theta0: "<<mylr.theta0<<std::endl;
            // std::cout<<"theta1: "<<mylr.theta1<<std::endl;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_delta = 0;
            for (uint32_t i = 0; i < length; i++) {
                T tmp_val;

                T pred = theta0 + theta1 * i;
                if (pred > ascii_vec_max[i]) {
                    tmp_val = pred - ascii_vec_max[i];
                    signvec.emplace_back(false);
                }
                else if (pred < ascii_vec_min[i]) {
                    tmp_val = ascii_vec_min[i] - pred;
                    signvec.emplace_back(true);
                }
                else {
                    tmp_val = 0;
                    signvec.emplace_back(false);
                }

                delta.emplace_back(tmp_val);

                if (tmp_val > max_delta) {
                    max_delta = tmp_val;
                }
            }

            uint32_t max_bit = 0;
            if (max_delta) {
                max_bit = bits_T(max_delta) + 1;
            }
            // std::cout<< "max_bit: " << max_bit << std::endl;
            memcpy(out, &max_bit, sizeof(uint32_t));
            out += sizeof(uint32_t);
            memcpy(out, &padding_length[nvalue], sizeof(uint8_t));
            out += sizeof(uint8_t);
            memcpy(out, &theta0, sizeof(T));
            out += sizeof(T);
            memcpy(out, &theta1, sizeof(T));
            out += sizeof(T);
            if (max_bit) {
                out = write_delta_string(delta.data(), signvec, string_wo_padding_length.data() + block_size * nvalue, out, max_bit, length);
            }

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


        bool pre_Bsearch(int left, int right, int* index, const char* key,
            bool* skip_linear_scan, int string_len, const char* firstkey_each_block,
            int search_len) {  // check which block the key belongs to
            while (left != right) {
                // The `mid` is computed by rounding up so it lands in (`left`, `right`].
                int64_t mid = left + (right - left + 1) / 2;


                // uint8_t first_key_length = reinterpret_cast<const uint8_t*>(firstkey_each_block+ mid * string_len)[0];
                // const size_t min_len = (first_key_length < search_len) ? first_key_length : search_len;
                // int cmp = memcmp(firstkey_each_block + mid * string_len + 1, key, min_len);
                // if (cmp == 0) {
                //   if (first_key_length < search_len)
                //     cmp = -1;
                //   else if (first_key_length > search_len)
                //     cmp = 1;
                // }
                int cmp = memcmp(firstkey_each_block + mid * string_len + 1, key, search_len);


                // std::string string_comp(firstkey_each_block + mid * string_len + 1, search_len);
                // std::cout<< search_len<<" "<<string_comp<<std::endl;
                // need another length info when storing each first key.
                if (cmp < 0) {
                    // Key at "mid" is smaller than "target". Therefore all
                    // blocks before "mid" are uninteresting.
                    left = mid;
                }
                else if (cmp > 0) {
                    // Key at "mid" is >= "target". Therefore all blocks at or
                    // after "mid" are uninteresting.
                    right = mid - 1;
                }
                else {

                    uint8_t first_key_length = reinterpret_cast<const uint8_t*>(firstkey_each_block + mid * string_len)[0];
                    if (first_key_length < search_len)
                        left = mid;
                    else if (first_key_length > search_len)
                        right = mid - 1;
                    else {
                        *skip_linear_scan = true;
                        left = right = mid;
                    }


                }
            }

            if (left == -1) {
                // All keys in the block were strictly greater than `target`. So the very
                // first key in the block is the final seek result.
                *skip_linear_scan = true;
                *index = 0;
            }
            else {
                *index = static_cast<uint32_t>(left);
            }
            return true;
        }

        // int Bsearch(int low, int high, int index, T& record)
        // {
        //     //std::cout<< "left: " << left << " right: " << right << std::endl;

        //     while (left != right) {
        //         // The `mid` is computed by rounding up so it lands in (`left`, `right`].
        //         int mid = left + (right - left + 1) / 2;
        //         int mid_block = mid / block_size_;

        //         int start_byte = reinterpret_cast<const int*>(data_ + data_offset + mid_block * 4)[0];

        //         T data_mid;
        //         uint8_t origin_string_length = randomdecodeArray8_string<T>(data_ + data_offset + block_number * 4 + start_byte, mid % block_size_, &data_mid, target.size());
        //         // std::string mid_string = convertToString<T>(&data_mid, origin_string_length, target.size());
        //         // std::cout<<"data_mid: "<<mid_string<<std::endl;
        //         if (data_mid == record) {

        //             if (origin_string_length == target.size()) {
        //                 left = right = mid;
        //                 *index = static_cast<uint32_t>(left);
        //                 return true;
        //             }
        //             else if (origin_string_length > target.size()) {
        //                 right = mid - 1;
        //             }
        //             else {
        //                 left = mid;
        //             }

        //         }
        //         else if (data_mid > record) {
        //             right = mid - 1;
        //         }
        //         else if (data_mid < record) {
        //             left = mid;
        //         }
        //     }
        //     *index = static_cast<uint32_t>(left) + 1;

        // }

        // bool TestBsearch(int sample_size, std::vector<std::string>& string_vec, int N)
        // {
        //     for (int i = 0; i < sample_size; i++)
        //     {
        //         int index = random(N);
        //         std::string tmpkey = string_vec[index];
        //         bool skip_linear = false;
        //         int index_search = 0;

        //         uint32_t N;
        //         N = reinterpret_cast<const uint32_t*>(data_)[0];
        //         int block_number = N / block_size_;
        //         while (block_size_ * block_number < N)
        //         {
        //             block_number++;
        //         }
        //         int data_offset = sizeof(uint32_t) * 4;
        //         uint8_t max_padding_length;

        //         int interval = (block_size_ / key_num_per_block_);
        //         int total_firstkey_num = N / interval;
        //         if (N % interval) {
        //             total_firstkey_num += 1;
        //         }
        //         if (padding_enable_) {
        //             max_padding_length = reinterpret_cast<const uint8_t*>(data_)[data_offset];
        //             data_offset += sizeof(uint8_t);
        //             data_offset += total_firstkey_num * (max_padding_length + sizeof(uint8_t));
        //         }
        //         if (padding_enable_) {
        //             pre_Bsearch(-1, total_firstkey_num - 1, &index_search, target.data(), &skip_linear,
        //                 max_padding_length + sizeof(uint8_t), data_ + sizeof(uint32_t) * 4 + sizeof(uint8_t), target.size());
        //         }
        //         int index_result = 0;
        //         if (!skip_linear) {
        //             T record = 0;
        //             convertToASCII_char<T>(target.data(), target.size(), &record);

        //             int left = index_search * interval;
        //             int block_left = left / block_size_;
        //             if (left % block_size_) {
        //                 block_left += 1;
        //             }
        //             int right = std::min(left + interval, (block_left + 1) * block_size_ - 1);
        //             if (index_search == total_firstkey_num - 1) {
        //                 right = N;
        //             }
        //             //std::cout<< "left: " << left << " right: " << right << std::endl;

        //             Bsearch(left, right, &index_result, record);

        //         }
        //         else {
        //             *index_result = index_search * interval;
        //         }

        //         assert(index_result == index);

        //     }
        //     return true;
        // }

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
