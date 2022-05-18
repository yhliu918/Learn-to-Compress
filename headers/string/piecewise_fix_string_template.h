
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

namespace Codecset
{
    template <typename T>
    class Leco_string
    {
    public:
        uint8_t *encodeArray8_string(std::vector<std::string> &string_vec, int start_idx, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint8_t *out = res;
            std::vector<T> ascii_vec;
            std::vector<long_int> long_int_vec; // because the calculation in linear regression may overflow, we choose to adopt long_int in the string_lr and convert its param to type T later
            std::vector<int> index;
            for (int i = 0; i < length; i++)
            {
                ascii_vec.emplace_back(convertToASCII<T>(string_vec[i + start_idx]));
                long_int_vec.emplace_back(convertToLongInt(string_vec[i + start_idx]));
                index.emplace_back(i);
            }

            string_lr mylr;
            mylr.caltheta(index, long_int_vec, length);

            T theta0 = 0;
            // auto theta0_len = (mpz_sizeinbase(mylr.theta0.backend().data(), 2) + 7) / 8;
            mpz_export(&theta0, NULL, -1, 1, 1, 0, mylr.theta0.backend().data());
            // print_u128_u(theta0);
            // printf("\n");

            T theta1 = 0;
            // auto theta1_len = (mpz_sizeinbase(mylr.theta1.backend().data(), 2) + 7) / 8;
            mpz_export(&theta1, NULL, -1, 1, 1, 0, mylr.theta1.backend().data());
            // print_u128_u(theta1);
            // printf("\n");
            // std::cout<<"theta0: "<<mylr.theta0<<std::endl;
            // std::cout<<"theta1: "<<mylr.theta1<<std::endl;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_delta = 0;
            for (auto i = 0; i < length; i++)
            {
                T tmp_val;
                if (ascii_vec[i] > theta1 * i + theta0)
                {
                    tmp_val = ascii_vec[i] - (theta1 * i + theta0);
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = (theta1 * i + theta0) - ascii_vec[i];
                    signvec.emplace_back(false); // means negative
                }

                delta.emplace_back(tmp_val);

                if (tmp_val > max_delta)
                {
                    max_delta = tmp_val;
                }
            }

            uint32_t max_bit = 0;
            if (max_delta)
            {
                max_bit = bits_T(max_delta) + 1;
            }
            // std::cout<< "max_bit: " << max_bit << std::endl;
            memcpy(out, &max_bit, sizeof(uint32_t));
            out += sizeof(uint32_t);

            memcpy(out, &theta0, sizeof(T));
            out += sizeof(T);
            memcpy(out, &theta1, sizeof(T));
            out += sizeof(T);
            if (max_bit)
            {
                out = write_delta_string(delta.data(), signvec, NULL,out, max_bit, length);
            }

            return out;
        }

        void decodeArray8_string(uint8_t *in, const size_t length, size_t nvalue, std::vector<std::string> &string_vec)
        {
            uint32_t maxbits;
            uint8_t *tmpin = in;
            memcpy(&maxbits, tmpin, sizeof(uint32_t));
            // std::cout<<"max bit "<<maxbits<<std::endl;
            tmpin += sizeof(uint32_t);

            T theta0;
            memcpy(&theta0, tmpin, sizeof(T));
            tmpin += sizeof(T);

            T theta1;
            memcpy(&theta1, tmpin, sizeof(T));
            tmpin += sizeof(T);

            if (maxbits == 0)
            {
                for (int i = 0; i < length; i++)
                {
                    T tmp_val = theta0 + theta1 * i;
                    string_vec.emplace_back(convertToString<T>(&tmp_val));
                }
            }
            else
            {
                read_all_fix_string<T>(tmpin, 0, 0, length, maxbits, theta1, theta0, string_vec);
            }
            return ;
        }

        void randomdecodeArray8(int block_ind, int idx, T *result)
        {

            uint8_t *tmpin = block_start_vec[block_ind];
            uint32_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint32_t));
            tmpin += sizeof(uint32_t);

            T theta0;
            memcpy(&theta0, tmpin, sizeof(T));
            tmpin += sizeof(T);

            T theta1;
            memcpy(&theta1, tmpin, sizeof(T));
            tmpin += sizeof(T);

            read_bit_fix_string<T>(tmpin, maxbits, idx, theta1, theta0, result, NULL);
            return;
        }

        bool pre_Bsearch(int left, int right, int *index, std::string &key, bool *skip_linear_scan, int string_len)
        { // check which block the key belongs to
            while (left != right)
            {
                // The `mid` is computed by rounding up so it lands in (`left`, `right`].
                int64_t mid = left + (right - left + 1) / 2;
                std::string tmp_key = firstkey_each_block.substr(string_len * mid, string_len);
                if (tmp_key < key)
                {
                    // Key at "mid" is smaller than "target". Therefore all
                    // blocks before "mid" are uninteresting.
                    left = mid;
                }
                else if (tmp_key > key)
                {
                    // Key at "mid" is >= "target". Therefore all blocks at or
                    // after "mid" are uninteresting.
                    right = mid - 1;
                }
                else
                {
                    *skip_linear_scan = true;
                    left = right = mid;
                }
            }

            if (left == -1)
            {
                // All keys in the block were strictly greater than `target`. So the very
                // first key in the block is the final seek result.
                *skip_linear_scan = true;
                *index = 0;
            }
            else
            {
                *index = static_cast<uint32_t>(left);
            }
            return true;
        }

        bool Bsearch(int low, int high, int index, T &key)
        {
            int mid;
            if (low > high)
            {
                return false;
            }
            mid = (low + high) / 2;
            T data_mid;
            randomdecodeArray8(index, mid % block_size, &data_mid);

            if (data_mid == key)
            {
                return true;
            }
            else if (data_mid > key)
            {
                return Bsearch(low, mid - 1, index, key);
            }
            else if (data_mid < key)
            {
                return Bsearch(mid + 1, high, index, key);
            }
            return false;
        }

        bool TestBsearch(int sample_size, std::vector<std::string> &string_vec, int N)
        {
            for (int i = 0; i < sample_size; i++)
            {
                int index = random(N);
                std::string tmpkey = string_vec[index];
                bool skip_linear = false;
                int index_search = 0;
                int string_len = tmpkey.size();
                int blocks = block_start_vec.size();
                pre_Bsearch(-1, blocks - 1, &index_search, tmpkey, &skip_linear, string_len);
                // std::cout<< "index_search "<<index_search<<std::endl;
                if (!skip_linear)
                {
                    T record = convertToASCII<T>(tmpkey);
                    int lower = index_search * block_size;
                    int upper = (index_search + 1) * block_size;
                    if (index_search == blocks - 1)
                    {
                        upper = N;
                    }
                    Bsearch(lower, upper, index_search, record);
                }

            }
            return true;
        }

        void setBlockSize(int block_size)
        {
            this->block_size = block_size;
            return;
        }

        void push_back_block(uint8_t *block)
        {
            block_start_vec.emplace_back(block);
            return;
        }

        void push_back_firstkey(std::string &firstkey)
        {
            firstkey_each_block.append(firstkey);
            return;
        }

    private:
        std::string firstkey_each_block;
        std::vector<uint8_t *> block_start_vec;
        int block_size;

    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
