
#ifndef PIECEWISEFIXOP_H_
#define PIECEWISEFIXOP_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {





    class piecewise_fix_op : public IntegerCODEC {
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
        // std::vector<int> mul_add_diff_set;
        // std::vector<int> mul_add_diff_set_minus;
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
            double theta1 = mylr.theta1;

            double pred = theta0;
            int max_error = 0;
            for (int i = 0;i < (long long)length;i++) {
                long long pred_mul = theta0 + theta1 * (double)i;
                long long pred_add = pred;
                if(pred_mul> pred_add){
                    mul_add_diff_set.push_back(i+nvalue*block_size);
                }
                if(pred_mul< pred_add){
                    mul_add_diff_set_minus.push_back(i+nvalue*block_size);
                }
                pred +=theta1;
                int tmp = (long long)in[i] - pred_mul;
                delta[i] = tmp;
                if (abs(tmp) > max_error) {
                    max_error = abs(tmp);
                }
            }

            int tmp_bit = 0;
            if (max_error) {
                tmp_bit = bits(max_error) + 1;
            }

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
            double theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            memcpy(&theta0, tmpin, sizeof(double));
            tmpin += sizeof(double);
            memcpy(&theta1, tmpin, sizeof(double));
            tmpin += sizeof(double);
            if (maxerror) {
                if (maxerror >= 31) {
                    read_all_default(tmpin, 0, 0, length, maxerror, theta1, theta0, out);
                }
                else {
                    read_all_bit_fix_add<uint32_t>(tmpin ,0,0, length, maxerror,theta1,theta0, out);
                    // read_all_bit_fix<uint32_t>(tmpin, 0, 0, length, maxerror, theta1, theta0, out);
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
            double theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            theta0 = reinterpret_cast<double*>(tmpin)[0];
            theta1 = reinterpret_cast<double*>(tmpin)[1];
            tmpin += sizeof(double) * 2;
            uint32_t tmp = 0;
            if (maxerror) {
                if (maxerror >= 31) {
                    //uint32_t * interpret = reinterpret_cast<uint32_t*>(tmpin);
                    tmp = read_bit_default(tmpin, maxerror, l, theta1, theta0, 0);
                    //return interpret[l];
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
            long long sum = 0;
            double theta0;
            double theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            memcpy(&theta0, tmpin, 8);
            tmpin += 8;
            memcpy(&theta1, tmpin, 8);
            tmpin += 8;
            double theta0_ = theta0;
            double theta1_ = theta1;
            // int64_t default_sum = 0;
            // for(int i=0;i<l;i++){
            //     // std::cout<<theta0_ + theta1_ * (double)i<< " "<<(long long)(theta0_ + theta1_ * (double)i)<<std::endl;
            //     default_sum += (long long)(theta0_ + theta1_ * (double)i);
            // }
            int64_t theta0_int = theta0;
            int64_t theta1_int = theta1;
            int64_t base_summation=0;
            // need to add counts of records that theta0 + i * theta1 <0
            // std::cout<<theta0_<<" "<<theta1_<<std::endl;
            if(theta0_<0 && theta1_ >0){
                // 0 ~ (-theta0_)/theta1_
                base_summation +=( (-theta0_)/theta1_);
            }
            if(theta1_ <0){
                //  (-theta0_)/theta1_ ~ l
                int add =  l - (-theta0_)/theta1_ ;
                if(add>0){
                    base_summation += add;
                }
                theta1_int --;
            }
            // std::cout<<"base_summation: "<<base_summation<<std::endl;
            base_summation += (theta0_int * (int)l + theta1_int * (int)((l - 1) * l / 2));
            
            theta0 -= theta0_int;
            theta1 -= theta1_int;


            std::vector<int> value_number_list;

            int origin_k = ((1.0 - theta0) / theta1);
            value_number_list.push_back(origin_k);

            int value = 1;
            int step = (1 / theta1);
            double epsilon = 1.0 - (double)step * theta1;
            int start_idx = origin_k + 1;
            int end_idx = 0;

            double residual = theta0 + theta1 + origin_k * theta1 - 1.0;

            double pred = theta0;
            int count = start_idx;

            while (count < l) {
                end_idx = start_idx + step;
                if (residual - epsilon >= 0) {
                    int number = end_idx - start_idx;
                    count += number;
                    if (count > l - 1) {
                        value_number_list.push_back(l - start_idx);
                        break;
                    }
                    else {
                        value_number_list.push_back(number);
                        residual -= epsilon;
                        value++;
                        start_idx = end_idx;
                    }
                }
                else {
                    int number = end_idx - start_idx + 1;
                    count += number;
                    if (count > l - 1) {
                        value_number_list.push_back(l - start_idx);
                        break;
                    }
                    else {
                        value_number_list.push_back(number);
                        residual = residual - epsilon + theta1;
                        value++;
                        start_idx = end_idx + 1;
                    }
                }
            }
       
            for (int i = 0;i < value_number_list.size();i++) {
                base_summation += value_number_list[i] * i;
            }
            // std::cout<<default_sum<<" "<<base_summation<<std::endl;
            // if(default_sum != base_summation){
            //     std::cout<<"error"<<std::endl;
            // }
            // assert(base_summation == default_sum);


            int64_t delta_sum = sum_all_deltas<uint32_t>(tmpin, l, maxerror);
            // int64_t ground_truth_sum = 0;
            // for(int i=0;i<l;i++){
            //     uint32_t decode = randomdecodeArray8(in, i, NULL, nvalue);
            //     int pred = (long long)(theta0_ + theta1_ * (double)i);
            //     int delta = decode - pred;
            //     ground_truth_sum+= delta;
            // }
           
            // assert(delta_sum == ground_truth_sum);
            
            delta_sum = delta_sum + base_summation;

            return delta_sum;
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
            return "piecewise_fix_op";
        }
        void destroy() {}
    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
