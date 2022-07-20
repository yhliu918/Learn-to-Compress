
#ifndef PIECEWISEFIX_STRING_SUBSET_TEMPLATE_H_
#define PIECEWISEFIX_STRING_SUBSET_TEMPLATE_H_

#include "../common.h"
#include "../codecs.h"
#include "../bit_read.h"
#include "../bit_write.h"
#include "lr_string.h"
#include "string_utils.h"
#include "bit_read_string.h"
#include "../headers/string/leco_uint256.h"
typedef leco_uint256 uint256_t;
#define INF 0x7f7fffff
#include <limits>

namespace Codecset
{
    inline int extract_common_prefix(std::vector<std::string>& string_vec, int start_ind,
        int block_length, int& common_prefix_length, std::string& common_prefix) {
        // compare the first and last string in each block, to extract common prefix
        std::string first_str = string_vec[start_ind];
        std::string last_str = string_vec[start_ind + block_length - 1];
        int compare_length = std::min(first_str.size(), last_str.size());
        common_prefix_length = 0;
        int j=0;
        for (j = 0; j < compare_length; j++) {
            if (first_str[j] != last_str[j]) {
                common_prefix_length = j;
                break;
            }
        }
        if(j==compare_length && common_prefix_length==0){
            common_prefix_length = compare_length;
        }
            
        common_prefix = first_str.substr(0, common_prefix_length);

        int max_padding_length = 0;
        for (int i = start_ind; i < start_ind + block_length; i++) {
            string_vec[i] = string_vec[i].substr(
                common_prefix_length, string_vec[i].size() - common_prefix_length);
            max_padding_length = std::max(max_padding_length, (int)string_vec[i].size());
        }
        return max_padding_length;  // block max key length(remain)
    }

    inline void Padd_one_string(std::string& s, char x, int padding_length) {
        s.append(padding_length - s.size(), x);
    }

    template <typename T>
    class Leco_string_subset {
    public:
        void Padding_string(std::vector<std::string>& string_vec, int start_ind,
            int block_length, char x, int padding_length, int min_ascii, int max_ascii) {
            totalsize_without_padding = 0;
            totalsize_with_padding = 0;
            for (int j = 0; j < block_length; j++) {
                totalsize_without_padding += string_vec[start_ind + j].size();
            }
            block_padding_length = padding_length;
            totalsize_with_padding += padding_length * block_length;
            for (int j = 0; j < block_length; j++) {
                int temp_len = string_vec[start_ind + j].size();
                string_wo_padding_length.emplace_back((uint8_t)temp_len);
                std::string tmp_str_max = string_vec[start_ind + j];
                std::string tmp_str_min = tmp_str_max;
                padding_max.emplace_back(tmp_str_max.append(padding_length - temp_len, max_ascii));
                padding_min.emplace_back(tmp_str_min.append(padding_length - temp_len, min_ascii));
                string_vec[start_ind + j].append(padding_length - temp_len, x);
            }
        }

        uint8_t* encodeArray8_string(std::vector<std::string>& string_vec,
            int start_idx, const size_t length, uint8_t* res,
            std::string* common_prefix, int common_prefix_length, uint8_t encoding_type) {

            uint8_t* out = res;
            // if(length == 1){
            //   out[0] = (uint8_t)string_vec[start_idx].size();
            //   out +=sizeof(uint8_t);
            //   T result = convertToASCII<T>(string_vec[start_idx])
            //   memcpy(out, &result, sizeof(T));
            //   out += sizeof(T);
            //   return out;
            // }                             

            std::vector<T> ascii_vec;
            std::vector<T> ascii_vec_min;
            std::vector<T> ascii_vec_max;
            std::vector<long_int>
                long_int_vec;  // because the calculation in linear regression may
            // overflow, we choose to adopt long_int in the string_lr
            // and convert its param to type T later
            std::vector<int> index;
            for (size_t i = 0; i < length; i++) {
                ascii_vec.emplace_back(convertToASCII_subset<T>(min_ascii, max_ascii, string_vec[i + start_idx]));
                long_int_vec.emplace_back(convertToLongInt_subset(min_ascii, max_ascii, string_vec[i + start_idx]));
                ascii_vec_min.emplace_back(convertToASCII_subset<T>(min_ascii, max_ascii, padding_min[i]));
                ascii_vec_max.emplace_back(convertToASCII_subset<T>(min_ascii, max_ascii, padding_max[i]));
                index.emplace_back(i);
            }

            string_lr mylr;
            mylr.caltheta(index, long_int_vec, length);

            T theta0 = 0;

            // auto theta0_len = (mpz_sizeinbase(mylr.theta0.backend().data(), 2) + 7) /
            // 8;
            mpz_export(&theta0, NULL, -1, 1, 1, 0, mylr.theta0.backend().data());
            

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
            int counter = 0;
            std::vector<uint8_t> string_wo_padding_length_vec;
            for (uint32_t i = 0; i < length; i++) {
                // when i is several times of interval, it is stored by first key (uncompressed format)
   
                T tmp_val;
                string_wo_padding_length_vec.emplace_back(string_wo_padding_length[i]);

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
                counter++;

                if (tmp_val > max_delta) {
                    max_delta = tmp_val;
                }
            }

            uint32_t max_bit = 0;
            if (max_delta) {
                max_bit = bits_T<T>(max_delta) + 1;
            }
            // std::cout<<"delta length "<<max_bit<<std::endl;

            memcpy(out, &encoding_type, sizeof(uint8_t));
            out += sizeof(uint8_t);
            memcpy(out, &common_prefix_length, sizeof(uint32_t));
            out += sizeof(uint32_t);
            memcpy(out, common_prefix->c_str(), common_prefix_length);
            // std::cout<<"common prefix "<<*common_prefix<<std::endl;
            out += common_prefix_length;

            // std::cout<< "max_bit: " << max_bit << std::endl;
            memcpy(out, &max_bit, sizeof(uint32_t));
            out += sizeof(uint32_t);
            memcpy(out, &block_padding_length, sizeof(uint8_t));
            out += sizeof(uint8_t);
            memcpy(out, &theta0, sizeof(T));
            out += sizeof(T);
            memcpy(out, &theta1, sizeof(T));
            out += sizeof(T);

            out = write_delta_string(
                delta.data(), signvec,
                string_wo_padding_length_vec.data(), out, max_bit,
                counter);


            return out;
        }

        void decodeArray8_string(uint8_t* in, const size_t length, size_t nvalue,
            std::vector<std::string>& string_vec) {
            uint32_t maxbits;
            uint8_t* tmpin = in;
            memcpy(&maxbits, tmpin, sizeof(uint32_t));
            // std::cout<<"max bit "<<maxbits<<std::endl;
            tmpin += sizeof(uint32_t);

            T theta0;
            memcpy(&theta0, tmpin, sizeof(T));
            tmpin += sizeof(T);

            T theta1;
            memcpy(&theta1, tmpin, sizeof(T));
            tmpin += sizeof(T);

            if (maxbits == 0) {
                for (int i = 0; i < length; i++) {
                    T tmp_val = theta0 + theta1 * i;
                    string_vec.emplace_back(convertToString_subset<T>(min_ascii, max_ascii, &tmp_val));
                }
            }
            else {
                read_all_fix_string<T>(tmpin, 0, 0, length, maxbits, theta1, theta0,
                    string_vec);
            }
            return;
        }



        void get_uncompressed_size(uint64_t& withpad, uint64_t& withoutpad) {
            withpad += totalsize_with_padding;
            withoutpad += totalsize_without_padding;
        }
        uint8_t get_max_padding_length() { return block_padding_length; }

        void write_string_wo_padding_len(std::string& buf) {
            for (char item : string_wo_padding_length) {
                buf.append(reinterpret_cast<const char*>(&item), sizeof(char));
            }
        }

        void write_block_padding_len(std::string& buf) {
            buf.append(reinterpret_cast<const char*>(&block_padding_length), sizeof(char));
        }

        void set_ascii_set(int min, int max){
            min_ascii = min;
            max_ascii = max;
        }

    private:
        uint64_t totalsize_without_padding;
        uint64_t totalsize_with_padding;
        std::vector<uint8_t> string_wo_padding_length;
        std::vector<std::string> padding_max;
        std::vector<std::string> padding_min;
        uint8_t block_padding_length;
        int min_ascii;
        int max_ascii;
    };

    template <typename T>
    std::string randomdecodeArray_leco_string(const uint8_t* in, int idx, T* result, int min_ascii, int max_ascii) {
        const uint8_t* tmpin = in;

        uint32_t maxbits = 0;
        memcpy(&maxbits, tmpin, sizeof(uint32_t));
        tmpin += sizeof(uint32_t);
        // std::cout<<"max bit "<<maxbits<<std::endl;
        uint8_t block_pad_length = 0;
        memcpy(&block_pad_length, tmpin, sizeof(uint8_t));
        tmpin += sizeof(uint8_t);

        T theta0 = 0;
        memcpy(&theta0, tmpin, sizeof(T));
        tmpin += sizeof(T);

        T theta1 = 0;
        memcpy(&theta1, tmpin, sizeof(T));
        tmpin += sizeof(T);

        uint8_t ori_length = 0;
        if(maxbits){
            read_bit_fix_string<T>(tmpin, maxbits, idx, theta1, theta0, result,
            &ori_length);
        }
        else{
            *result = theta0 + theta1 * idx;
            // ori_length = block_pad_length;
        }
        std::string tmp = convertToString_subset<T>(min_ascii, max_ascii, result, block_pad_length);
        tmp = tmp.substr(0, ori_length);
        // *result = *result >> (uint8_t)(8 * shift);
        // *result = *result << (uint8_t)(8 * shift);

        return tmp;
    }


    std::string randomdecode_string(uint8_t* in, int idx, int min_ascii, int max_ascii) {
        uint8_t* tmpin = in;

        uint8_t encoding_type = tmpin[0];
        tmpin++;
        int common_prefix_length = reinterpret_cast<const int*>(tmpin)[0];
        tmpin += sizeof(int);

        std::string common_prefix(reinterpret_cast<char*>(tmpin), common_prefix_length);
        tmpin += common_prefix_length;
        
        switch(encoding_type){
            case 0:{
                uint64_t result;
                std::string suffix =randomdecodeArray_leco_string<uint64_t>(tmpin, idx, &result, min_ascii, max_ascii);
                common_prefix.append(suffix);
                // common_prefix.append(convertToString_subset<uint64_t>(min_ascii, max_ascii, &result, ori_length));
                return common_prefix;
            }
            case 1: {
                uint128_t result;
                std::string suffix =randomdecodeArray_leco_string<uint128_t>(tmpin, idx, &result, min_ascii, max_ascii);
                common_prefix.append(suffix);
                // common_prefix.append(convertToString_subset<uint128_t>(min_ascii, max_ascii, &result, ori_length));
                return common_prefix;
            }
            case 2:{
                uint256_t result;
                std::string suffix =randomdecodeArray_leco_string<uint256_t>(tmpin, idx, &result, min_ascii, max_ascii);
                common_prefix.append(suffix);
                // common_prefix.append(convertToString_subset<uint256_t>(min_ascii, max_ascii, &result, ori_length));
                return common_prefix;
            }
            default:
                std::cout<<"Mpz_t not implemented yet, unknown encoding type"<<std::endl;
                return "";


        }
        
    }
    
} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
