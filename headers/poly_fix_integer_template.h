
#ifndef POLY_INTEGER_TEMPLATE_H_
#define POLY_INTEGER_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "bit_write.h"
#include "lr.h"
#include "bit_read.h"
#include "polynomial_regression.hpp"
#define INF 0x7f7fffff

namespace Codecset
{
    template <typename T, int degree>
    class Leco_int_poly
    {
    public:
        std::vector<std::pair<int, int>> mul_add_diff_set;
        int blocks;
        int block_size;

        void init( int block, int block_s)
        {
            blocks = block;
            block_size = block_s;
        }
        uint8_t *encodeArray8_int(const T *data, const size_t length, uint8_t *res, size_t nvalue)
        {
            uint8_t *out = res;
            std::vector<T> data_vec = std::vector<T>(data, data + length);
            auto simple_fixed = andviane::polynomial_regression<degree>(data_vec); 
            
            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;
            for (auto i = 0; i < length; i++)
            {
                T tmp_val;
                int64_t pred = simple_fixed(i);

                if ((int64_t)data[i] > pred)
                {
                    tmp_val = (int64_t)data[i] - pred;
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = pred - (int64_t)data[i];
                    signvec.emplace_back(false); // means negative
                }
                delta.emplace_back(tmp_val);

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }
            uint8_t max_bit = 0;
            if (max_error)
            {
                max_bit = bits_int_T(max_error) + 1;
            }

            if (max_bit > sizeof(T) * 8)
            {
                max_bit = sizeof(T) * 8;
            }
            memcpy(out, &max_bit, sizeof(max_bit));
            out += sizeof(max_bit);
            if (max_bit == sizeof(T) * 8)
            {
                for (auto i = 0; i < length; i++)
                {
                    memcpy(out, &data[i], sizeof(T));
                    out += sizeof(T);
                }
                return out;
            }
            // memcpy(out, &degree, sizeof(uint8_t));
            for(auto a: simple_fixed){
                memcpy(out, &a, sizeof(double));
                out += sizeof(double);
            }
            if (max_bit)
            {
                out = write_delta_int_T(delta, signvec, out, max_bit, length);
            }

            return out;
        }

        T *decodeArray8(const uint8_t *in, const size_t length, T *out, size_t nvalue)
        {
            T *res = out;
            // start_index + bit + theta0 + theta1 + numbers + delta

            const uint8_t *tmpin = in;

            uint8_t maxerror = tmpin[0];
            tmpin++;
            if (maxerror == 127)
            {
                T tmp_val;
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                return out;
            }
            if (maxerror == 126)
            {
                T tmp_val;
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                memcpy(&tmp_val, tmpin + sizeof(T), sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                return out;
            }
            if (maxerror >= sizeof(T) * 8 - 1)
            {
                // out = reinterpret_cast<T*>(tmpin);
                memcpy(out, tmpin, sizeof(T) * length);
                return out;
                // read_all_default(tmpin, 0, 0, length, maxerror, theta1, theta0, res);
            }
            // int degree = 0;
            // memcpy(&degree, tmpin, sizeof(uint8_t));
            std::array<double, degree + 1> theta;
            for (auto i = 0; i < degree + 1; i++)
            {
                memcpy(&theta[i], tmpin, sizeof(double));
                tmpin += sizeof(double);
            }
            andviane::Polynomial<degree> poly(theta);

            if (maxerror)
            {

                // read_all_bit_fix_add<T>(tmpin, 0, 0, length, maxerror, theta1, theta0, res);
                read_all_bit_fix_poly<T,degree,double>(tmpin, 0, 0, length, maxerror, res, &poly);
            }
            else
            {
                for (int j = 0; j < length; j++)
                {
                    res[j] = (long long)(poly(j));
                }
                // double pred = theta0;
                // for (int i = 0;i < length;i++) {
                //     res[i] = (long long)pred;
                //     pred += theta1;
                // }
            }

            return out;
        }

        int filter_range(uint8_t *in, const size_t length, T filter, uint32_t *out, int block_id)
        {
            // only consider > filter, return [return_value, end] for demo
            int block_start = block_id * block_size;
            int counter = 0;
            double theta0;
            double theta1;
            uint8_t maxerror;
            uint8_t *tmpin = in;
            maxerror = tmpin[0];
            tmpin++;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);
            int64_t delta_interval = 0;
            if (maxerror)
            {
                delta_interval = (1L << (maxerror - 1));
            }
            int thre = std::max((double)(filter + 1 - delta_interval - theta0) / theta1, 0.);
            if (thre >= length)
            {
                return counter;
            }
            else
            {
                if (maxerror)
                {
                    counter = read_all_bit_fix_range<T>(tmpin, 0, thre, length, maxerror, theta1, theta0, out, filter, block_start);
                }
                else
                {
                    for (int i = thre; i < length; i++)
                    {
                        long long pred = (long long)(theta0 + theta1 * (double)i);
                        if (pred > filter)
                        {
                            out[0] = block_start + i;
                            out++;
                            counter++;
                        }
                    }
                }
            }
            return counter;
        }

        int filter_range_close(uint8_t *in, const size_t length, uint32_t *out, int block_id, T filter1, T filter2, int base)
        {
            // only  filter2 > consider > filter1, return [return_value, end] for demo
            // (base*i + filter1, base*i + filter2]
            int block_start = block_id * block_size;
            int counter = 0;
            double theta0;
            double theta1;
            uint8_t maxerror;
            uint8_t *tmpin = in;
            maxerror = tmpin[0];
            tmpin++;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);
            int64_t delta_interval = 0;
            if (maxerror)
            {
                delta_interval = (1L << (maxerror - 1));
            }
            int thre1 = std::max((double)(filter1 - delta_interval - theta0) / theta1 -1, 0.);
            int thre2 = std::min((double)(filter2 + delta_interval - theta0) / theta1 +1, (double)length);
            while(thre1 <length && thre2 > 0){
                if (maxerror)
                {
                    counter += read_all_bit_fix_range_close<T>(tmpin, 0, thre1, thre2, length, maxerror, theta1, theta0, out, filter1, filter2, block_start);
                }
                else
                {
                    for (int i = thre1; i <= thre2; i++)
                    {
                        long long pred = (long long)(theta0 + theta1 * (double)i);
                        if (pred > filter1 && pred < filter2)
                        {
                            out[0] = block_start + i;
                            out++;
                            counter++;
                        }
                    }
                }
                filter1 += base;
                filter2 += base;
                thre1 = std::max((double)(filter1 - delta_interval - theta0) / theta1 -1, 0.);
                thre2 = std::min((double)(filter2 + delta_interval - theta0) / theta1 +1, (double)length);

            }
          
            return counter;
        }

        T randomdecodeArray8(const uint8_t *in, int to_find, uint32_t *out, size_t nvalue)
        {

            const uint8_t *tmpin = in;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);

            if (maxbits == sizeof(T) * 8)
            {
                T tmp_val = reinterpret_cast<const T *>(tmpin)[to_find];
                return tmp_val;
            }
            // int degree = 0;
            // memcpy(&degree, tmpin, sizeof(uint8_t));
            std::array<double, degree + 1> theta;
            for (auto i = 0; i < degree + 1; i++)
            {
                memcpy(&theta[i], tmpin, sizeof(double));
                tmpin += sizeof(double);
            }
            andviane::Polynomial<degree> poly(theta);
            T tmp_val;
            if (maxbits)
            {
                tmp_val = read_bit_fix_int_wo_round_poly<T,degree,double>(tmpin, maxbits, to_find, &poly);
            }
            else
            {
                tmp_val = poly(to_find);
            }
            return tmp_val;
        }
    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
