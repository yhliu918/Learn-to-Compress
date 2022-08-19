
#ifndef PIECEWISEFIX_INTEGER_TEMPLATE_H_
#define PIECEWISEFIX_INTEGER_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "bit_write.h"
#include "lr.h"
#include "bit_read.h"
#define INF 0x7f7fffff

namespace Codecset
{
    template <typename T>
    class Leco_int
    {
    public:
        std::vector<std::pair<int,int>> mul_add_diff_set;
        int blocks;
        int block_size;
        
        void init(int block, int block_s){
            blocks = block;
            block_size = block_s;
        }
        uint8_t* encodeArray8_int(T* data, const size_t length, uint8_t* res, size_t nvalue)
        {
            uint8_t* out = res;
            std::vector<double> indexes;
            std::vector<double> keys;
            for (uint32_t i = 0; i < length; i++) {
                indexes.emplace_back((double)i);
                keys.emplace_back((double)data[i]);
            }

            lr mylr;
            mylr.caltheta(indexes, keys, length);
            double theta0 = mylr.theta0;
            double theta1 = mylr.theta1;

            // int128_t max_error_delta = INT64_MIN;
            // int128_t min_error_delta = INT64_MAX;
            // for (auto i = 0; i < length; i++){
            //     int128_t tmp_val =  (int128_t)data[i] - (int128_t)(theta0 + theta1 * (double)i);
            //     if(tmp_val>max_error_delta) max_error_delta = tmp_val;
            //     if(tmp_val<min_error_delta) min_error_delta = tmp_val;
            // }
            // theta0+=(max_error_delta+min_error_delta)/2.;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;
            for (auto i = 0; i < length; i++)
            {
                T tmp_val;
                if (data[i] > (T)(theta0 + theta1 * (double)i))
                {
                    tmp_val = data[i] - (T)(theta0 + theta1 * (double)i);
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = (T)(theta0 + theta1 * (double)i) - data[i];
                    signvec.emplace_back(false); // means negative
                }

                delta.emplace_back(tmp_val);

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }

            // double pred = theta0;
            // T max_error = 0;
            // for (auto i = 0; i < length; i++)
            // {
            //     T pred_mul = (theta0 + theta1 * (double)i);
            //     T pred_add = pred;
            //     if(pred_mul> pred_add){
            //         mul_add_diff_set.push_back(std::make_pair(i+nvalue*block_size,pred_mul - pred_add));
            //     }
            //     if(pred_mul< pred_add){
            //         mul_add_diff_set.push_back(std::make_pair(i+nvalue*block_size,-(int)(pred_add - pred_mul)));
            //     }
            //     pred +=theta1;

            //     T tmp_val;
            //     if (data[i] > pred_mul)
            //     {
            //         tmp_val = data[i] - pred_mul;
            //         signvec.emplace_back(true); // means positive
            //     }
            //     else
            //     {
            //         tmp_val = pred_mul - data[i];
            //         signvec.emplace_back(false); // means negative
            //     }

            //     delta.emplace_back(tmp_val);

            //     if (tmp_val > max_error)
            //     {
            //         max_error = tmp_val;
            //     }
            // }

            uint8_t max_bit = 0;
            if (max_error)
            {
                max_bit = bits_int_T(max_error) + 1;
            }
            // std::cout<< "max_bit: " << (int)max_bit << std::endl;
            if (max_bit > sizeof(T) * 8) {
                max_bit = sizeof(T) * 8;
            }
            memcpy(out, &max_bit, sizeof(max_bit));
            out += sizeof(max_bit);
            if (max_bit == sizeof(T) * 8) {
                for (auto i = 0; i < length; i++)
                {
                    memcpy(out, &data[i], sizeof(T));
                    out += sizeof(T);
                }
                return out;
            }

            memcpy(out, &theta0, sizeof(double));
            out += sizeof(double);
            memcpy(out, &theta1, sizeof(double));
            out += sizeof(double);
            if (max_bit)
            {
                out = write_delta_int_T(delta, signvec, out, max_bit, length);
            }

            return out;
        }


        T* decodeArray8(uint8_t *in, const size_t length, T* out, size_t nvalue) {
            T* res = out;
            //start_index + bit + theta0 + theta1 + numbers + delta
            double theta0;
            double theta1;
            uint8_t* tmpin = in;

            uint8_t maxerror = tmpin[0];
            tmpin++;
            if (maxerror == 127) {
                T tmp_val;
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                return out;
            }
            if (maxerror == 126) {
                T tmp_val;
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                memcpy(&tmp_val, tmpin + sizeof(T), sizeof(tmp_val));
                res[0] = tmp_val;
                res++;
                return out;
            }


            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);

            if (maxerror) {
                if (maxerror >= sizeof(T) * 8 - 1) {
                    // read_all_default(tmpin, 0, 0, length, maxerror, theta1, theta0, res);
                }
                else {
                    // read_all_bit_fix_add<T>(tmpin, 0, 0, length, maxerror, theta1, theta0, res);
                    read_all_bit_fix<T>(tmpin, 0, 0, length, maxerror, theta1, theta0, res);
                }
            }
            else {
                for (int j = 0;j < length;j++) {
                    res[j] = (long long)(theta0 + theta1 * (double)j);
                }
                // double pred = theta0;
                // for (int i = 0;i < length;i++) {
                //     res[i] = (long long)pred;
                //     pred += theta1;
                // }
            }

            return out;
        }



        T randomdecodeArray8(uint8_t* in, int to_find, uint32_t* out, size_t nvalue)
        {

            uint8_t* tmpin = in;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);

            if (maxbits == sizeof(T) * 8) {
                T tmp_val = reinterpret_cast<T*>(tmpin)[to_find];
                return tmp_val;
            }

            double theta0;
            memcpy(&theta0, tmpin, sizeof(double));
            tmpin += sizeof(double);

            double theta1;
            memcpy(&theta1, tmpin, sizeof(double));
            tmpin += sizeof(double);


            T tmp_val;
            if (maxbits) {
                tmp_val = read_bit_fix_int_wo_round<T>(tmpin, maxbits, to_find, theta1, theta0);
            }
            else {
                tmp_val = (theta0 + (double)to_find * theta1);
            }
            return tmp_val;
        }



    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
