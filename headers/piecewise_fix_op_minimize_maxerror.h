
#ifndef PIECEWISEFIXOP_MAX_H_
#define PIECEWISEFIXOP_MAX_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {





    class piecewise_fix_op_max : public IntegerCODEC {
    public:
        using IntegerCODEC::encodeArray;
        using IntegerCODEC::decodeArray;
        using IntegerCODEC::randomdecodeArray;
        using IntegerCODEC::encodeArray8;
        using IntegerCODEC::decodeArray8;
        using IntegerCODEC::randomdecodeArray8;
        using IntegerCODEC::init;
        using IntegerCODEC::summation;


        int block_num;
        int block_size;
        //double * models;

        void init(int blocks, int blocksize, int extra) {
            block_num = blocks;
            block_size = blocksize;
            //models = new double [block_num*2];
        }

        int random(int m) {
            return rand() % m;
        }

        // bit + theta0 + theta1 + delta   
        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* res, size_t nvalue) {

            uint8_t* out = res;
            int* delta = new int[length];

            lr_int mylr;
            mylr.caltheta(in, length);


            double theta0 = mylr.theta0;
            float theta1 = mylr.theta1;


            int max_error_delta = INT_MIN;
            int min_error_delta = INT_MAX;
            for (int i = 0;i < (long long)length;i++) {
                int tmp = in[i] - (long long)(theta0 + theta1 * (double)i);
                if (tmp > max_error_delta) {
                    max_error_delta = tmp;
                }
                if (tmp < min_error_delta) {
                    min_error_delta = tmp;
                }
            }
            theta0 += (max_error_delta + min_error_delta) / 2.0;


            double pred = theta0;
            int max_error = 0;
            for (int i = 0;i < (long long)length;i++) {
                long long pred_mul = theta0 + theta1 * (double)i;
                long long pred_add = pred;
                if (pred_mul > pred_add) {
                    mul_add_diff_set.push_back(i + nvalue * block_size);
                }
                if (pred_mul < pred_add) {
                    mul_add_diff_set_minus.push_back(i + nvalue * block_size);
                }
                pred += theta1;
                int tmp = (long long)in[i] - pred_mul;
                delta[i] = tmp;
                if (abs(tmp) > max_error) {
                    max_error = abs(tmp);
                }
            }
            // std::cout<<nvalue*block_size<<" "<<nvalue*block_size+length-1<<" "<<theta0<<" "<<theta1<<std::endl;



            int tmp_bit = 0;
            if (max_error) {
                tmp_bit = bits(max_error) + 1;
            }

            //std::cout<<"bit_length: "<<tmp_bit<<std::endl;
            out[0] = (uint8_t)tmp_bit;
            out++;

            memcpy(out, &theta0, sizeof(theta0));
            out += sizeof(theta0);
            memcpy(out, &theta1, sizeof(theta1));
            out += sizeof(theta1);

            if (tmp_bit) {
                if (tmp_bit >= 31) {
                    out = write_delta_default(in, out, 32, length);
                }
                else {
                    out = write_delta_T(delta, out, tmp_bit, length);
                }
            }

            free(delta);

            return out;

        }


        uint32_t* decodeArray8(uint8_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            double theta0;
            float theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            maxerror = tmpin[0];
            tmpin++;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);
            if (maxerror) {
                if (maxerror >= 31) {
                    read_all_default(tmpin, 0, 0, length, maxerror, theta1, theta0, out);
                }
                else {
                    // read_all_bit_fix<uint32_t>(tmpin, 0, 0, length, maxerror, theta1, theta0, out);
                    read_all_bit_fix_add<uint32_t>(tmpin ,0,0, length, maxerror,theta1,theta0, out);
                }
            }
            else {
                // for (int i = 0;i < length;i++) {
                //     out[i] = (long long)(theta0 + theta1 * (double)i);
                // }
                double pred = theta0;
                for (int i = 0;i < length;i++) {
                    out[i] = (long long)pred;
                    pred += theta1;
                }
            }


            return out;
        }

        uint32_t randomdecodeArray8(uint8_t* in, const size_t l, uint32_t* out, size_t nvalue) {
            double theta0;
            float theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            theta0 = reinterpret_cast<double*>(tmpin)[0];
            tmpin += sizeof(theta0);
            theta1 = reinterpret_cast<float*>(tmpin)[0];
            tmpin += sizeof(theta1);

            uint32_t tmp = 0;
            if (maxerror) {
                if (maxerror >= 31) {
                    tmp = read_bit_default(tmpin, maxerror, l, theta1, theta0, 0);
                }
                else {
                    tmp = read_bit_fix_T(tmpin, maxerror, l, theta1, theta0, 0);
                }
            }
            else {
                tmp = (long long)(theta0 + theta1 * (double)l);
            }

            return tmp;

        }
        uint64_t summation(uint8_t* in, const size_t l, size_t nvalue) {
            return 0;
        }

        uint32_t* encodeArray(uint32_t* in, const size_t length, uint32_t* out,
            size_t nvalue) {
            std::cout << "Haven't implement. Please try uint8_t one..." << std::endl;
            return out;
        }
        uint32_t* decodeArray(uint32_t* in, const size_t length,
            uint32_t* out, size_t nvalue) {
            std::cout << "Haven't implement. Please try uint8_t one..." << std::endl;
            return out;
        }
        uint32_t randomdecodeArray(uint32_t* in, const size_t l, uint32_t* out, size_t nvalue) {
            std::cout << "Haven't implement. Please try uint8_t one..." << std::endl;
            return 1;
        }
        uint32_t get_block_nums() {
            return 1;
        }
        std::string name() const {
            return "piecewise_fix_op_max";
        }
        void destroy() {}
    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
