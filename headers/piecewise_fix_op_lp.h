
#ifndef PIECEWISEFIXOP_LP_H_
#define PIECEWISEFIXOP_LP_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {

    class piecewise_fix_op_lp : public IntegerCODEC {
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
        std::string filename = "newman.log";
        std::vector<float> theta0_vec;
        std::vector<float> theta1_vec;
        //double * models;

        void init(int blocks, int blocksize, int extra) {
            block_num = blocks;
            block_size = blocksize;
            std::ifstream srcFile("/home/lyh/Learn-to-Compress/scripts/" + filename, std::ios::in);
            while (srcFile.good())
            {
                float tmp0;
                float tmp1;
                
                srcFile >> tmp0 >> tmp1;
                theta0_vec.push_back(tmp0);
                theta1_vec.push_back(tmp1);   
                if (!srcFile.good())
                {
                    break;
                }
            }
            srcFile.close();
            //models = new double [block_num*2];
        }

        int random(int m) {
            return rand() % m;
        }

        // bit + theta0 + theta1 + delta   
        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* res, size_t nvalue) {
            uint8_t* out = res;
            int* delta = new int[length];
            int max_error = 0;
            float theta0 = theta0_vec[nvalue];
            float theta1 = theta1_vec[nvalue];
            if(theta1<0.00001)
            {
                theta1 = 0.0;
            }
            for (int i = 0;i < (long long)length;i++) {
                int tmp = (long long)in[i] - (long long)(theta0 + theta1 * (double)i);
                delta[i] = tmp;
                if (abs(tmp) > max_error) {
                    max_error = abs(tmp);
                }
            }

            int tmp_bit = 0;
            if (max_error) {
                tmp_bit = bits(max_error) + 1;
            }
            // std::cout<<tmp_bit<<std::endl;
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
            float theta0;
            float theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            memcpy(&maxerror, tmpin, 1);
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

                    read_all_bit_fix<uint32_t>(tmpin, 0, 0, length, maxerror, theta1, theta0, out);
                }
            }
            else {
                for (int i = 0;i < length;i++) {
                    out[i] = (long long)(theta0 + theta1 * (double)i);
                }
            }


            return out;
        }

        uint32_t randomdecodeArray8(uint8_t* in, const size_t l, uint32_t* out, size_t nvalue) {
            float theta0;
            float theta1;
            uint8_t maxerror;
            uint8_t* tmpin = in;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);
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
        /*
        uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
            long long sum =0;
            double theta0;
            double theta1;
            uint8_t maxerror;
            uint8_t * tmpin=in;
            memcpy(&maxerror,tmpin,1);
            tmpin++;
            memcpy(&theta0,tmpin,8);
            tmpin+=8;
            memcpy(&theta1,tmpin,8);
            tmpin+=8;
            if(maxerror>=31){
                sum = sum_all_default(tmpin ,0,0, nvalue, maxerror);
            }
            else{
                //uint32_t total = nvalue *(nvalue-1)/2;
                sum = sum_all_bit_fix(tmpin ,0,0, nvalue, maxerror,theta1);
                sum +=((long long) theta0 )*nvalue;
                //sum +=(long long)(theta1*total);
            }
            return (uint64_t)sum;
        }
        */
       uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
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
            return "piecewise_fix_op_lp";
        }
        void destroy() {}
    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
