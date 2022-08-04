
#ifndef PIECEWISE_FITING_H_
#define PIECEWISE_FITING_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {

    class piecewise_fiting : public IntegerCODEC {
    public:
        using IntegerCODEC::encodeArray;
        using IntegerCODEC::decodeArray;
        using IntegerCODEC::randomdecodeArray;
        using IntegerCODEC::encodeArray8;
        using IntegerCODEC::decodeArray8;
        using IntegerCODEC::randomdecodeArray8;
        using IntegerCODEC::init;

        std::vector<uint8_t*> block_start_vec;
        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        uint64_t total_byte = 0;
        int overhead = 10;
        uint32_t* array;
        // int tolerance = 0;
        int block_num;
        int block_size;

        //start_index + bit + theta0 + theta1 + numbers + delta
        void init(int blocks, int blocksize, int delta) {
            block_num = blocks;
            block_size = blocksize;
            overhead = delta; // add some punishing item

        }
        uint32_t lower_bound(uint32_t v, uint32_t len)
        {
            uint32_t m;
            uint32_t x = 0;
            uint32_t y = len - 1;
            while (x <= y)
            {

                m = x + (y - x) / 2;
                if (v < segment_index[m]) y = m - 1;
                else x = m + 1;
            }
            return y;

        }

        void newsegment(uint32_t origin_index, uint32_t end_index) {
            if(end_index - origin_index + 1 ==2){
                return newsegment_2(origin_index,end_index);
            }
            if(end_index - origin_index + 1 ==1){
                return newsegment_1(origin_index,end_index);
            }
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t) * 2);
            uint8_t* out = descriptor;
            int length = end_index - origin_index + 1;
            std::vector<double> indexes;
            std::vector<double> keys;
            for (int j = origin_index;j <= end_index;j++) {
                indexes.emplace_back(j - origin_index);
                keys.emplace_back(array[j]);
            }

            lr mylr;
            mylr.caltheta(indexes, keys, length);
            double final_slope = mylr.theta1;
            double theta0 = mylr.theta0;

            uint32_t final_max_error = 0;
            int* delta_final = new int[end_index - origin_index + 1];

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                long long  pred = (long long)(theta0 + final_slope * (double)(j - origin_index));

                uint32_t tmp_error = abs(pred - array[j]);
                delta_final[j - origin_index] = array[j] - pred;
                if (tmp_error > final_max_error) {
                    final_max_error = tmp_error;
                }
            }
            uint32_t delta_final_max_bit = 0;
            if (final_max_error) {
                delta_final_max_bit = bits_int_T<uint32_t>(final_max_error) + 1;
            }



            if (delta_final_max_bit >= 32) {
                delta_final_max_bit = 32;
            }

            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)delta_final_max_bit;
            out++;

            memcpy(out, &theta0, sizeof(theta0));
            out += sizeof(theta0);

            memcpy(out, &final_slope, sizeof(final_slope));
            out += sizeof(final_slope);
            if (delta_final_max_bit) {
                if (delta_final_max_bit == 32) {
                    out = write_delta_default(array + origin_index, out, delta_final_max_bit, end_index - origin_index + 1);
                }
                else {
                    out = write_delta_T(delta_final, out, delta_final_max_bit, (end_index - origin_index + 1));

                }
            }



            delete[] delta_final;


            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);
            segment_length.push_back(segment_size);
            total_byte += segment_size;
            // if(origin_index == 2024){
            //     std::cout<<segment_size<<" "<<end_index<<std::endl;
            // }
        }

        void newsegment_2(uint32_t origin_index, uint32_t end_index) {
            // if(origin_index==1636 || origin_index+1 == 1636){
            //     std::cout<<"hello"<<std::endl;
            // }
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t));
            uint8_t* out = descriptor;
            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)126; // this means that this segment only has two points
            out++;
            memcpy(out, &array[origin_index], sizeof(uint32_t));
            out += sizeof(uint32_t);
            memcpy(out, &(array[origin_index + 1]), sizeof(uint32_t));
            out += sizeof(uint32_t);

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);
            segment_length.push_back(segment_size);

            total_byte += segment_size;
        }

        void newsegment_1(uint32_t origin_index, uint32_t end_index) {
            // if(origin_index == 1636){
            //     std::cout<<origin_index<<std::endl;
            // }
            uint8_t* descriptor = (uint8_t*)malloc(3 * sizeof(uint64_t));
            uint8_t* out = descriptor;
            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)127; // this means that this segment only has one point
            out++;
            memcpy(out, &array[origin_index], sizeof(uint32_t));
            out += sizeof(uint32_t);

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_length.push_back(segment_size);
            segment_index.push_back(origin_index);

            total_byte += segment_size;
        }


        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* res, size_t nvalue) {
            array = in;
            std::vector<uint32_t> indexes;
            for (uint32_t i = 0; i < nvalue; i++) {
                indexes.push_back(i);
            }
            double high_slope = (double)INF;
            double low_slope = 0.;
            long long origin_key = in[0];
            int origin_index = indexes[0];
            int end_index = indexes[0];
            int tmp_delta_bit = 0;
            int tmp_max_delta = 0;
            for (int i = 1; i < (long long)nvalue; i++) {
                long long key = in[i];
                int id = indexes[i];
                if (i == nvalue - 1) {
                    newsegment(origin_index, id);
                    break;
                }
                double tmp_point_slope = ((key - origin_key) + 0.0) / ((id - origin_index) + 0.0);
                if (tmp_point_slope >= low_slope && tmp_point_slope <= high_slope) {
                    double tmp_high_slope = ((key + overhead - origin_key) + 0.0) / ((id - origin_index) + 0.0);
                    double tmp_low_slope = ((key - overhead - origin_key) + 0.0) / ((id - origin_index) + 0.0);
                    if (tmp_high_slope <= high_slope) {
                        high_slope = tmp_high_slope;
                    }
                    if (low_slope <= tmp_low_slope) {
                        low_slope = tmp_low_slope;
                    }
                    end_index = id;
                }
                else {
                    int max_error = 0;
                    newsegment(origin_index, end_index);

                    
                    high_slope = (double)INF;
                    low_slope = 0.0;
                    origin_index = id;
                    origin_key = key;
                    end_index = id;


                }

            }

            return res;
        }
        

        uint32_t* decodeArray8(uint8_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            uint32_t* res = out;
            //start_index + bit + theta0 + theta1 + numbers + delta
            segment_index.push_back(length);
            double theta0;
            double theta1;
            uint8_t maxerror;
            for (int i = 0;i < block_start_vec.size();i++) {
                int segment_length = segment_index[i + 1] - segment_index[i];
                uint8_t* tmpin = block_start_vec[i];
                tmpin += sizeof(uint32_t);
                maxerror = tmpin[0];
                tmpin++;
                if (maxerror == 127) {
                    uint32_t tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }
                if (maxerror == 126) {
                    uint32_t tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    memcpy(&tmp_val, tmpin + sizeof(uint32_t), sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }


                memcpy(&theta0, tmpin, sizeof(theta0));
                tmpin += sizeof(theta0);
                memcpy(&theta1, tmpin, sizeof(theta1));
                tmpin += sizeof(theta1);
                if (maxerror) {
                    if (maxerror >= sizeof(uint32_t) * 8 - 1) {
                        read_all_default(tmpin, 0, 0, segment_length, maxerror, theta1, theta0, res);
                    }
                    else {
                        read_all_bit_fix<uint32_t>(tmpin, 0, 0, segment_length, maxerror, theta1, theta0, res);
                    }
                }
                else {
                    for (int j = 0;j < segment_length;j++) {
                        res[j] = (long long)(theta0 + theta1 * (double)j);
                    }
                }
                res += segment_length;
            }
            return out;
        }



    uint32_t randomdecodeArray8(uint8_t* in, const size_t l, uint32_t* out, size_t nvalue) {

        uint32_t length = segment_index.size();
        uint8_t* this_block = block_start_vec[lower_bound(l, length)];

        uint8_t* tmpin = this_block;
        double theta0;
        double theta1;
        uint8_t maxerror;
        uint32_t start_ind;
        uint32_t tmp = 0;
        memcpy(&start_ind, tmpin, 4);
        tmpin += 4;
        memcpy(&maxerror, tmpin, 1);
        tmpin++;


        if (maxerror == 127) {
            memcpy(&tmp, tmpin, 4);
            return tmp;
        }
        if (maxerror == 126) {
            if (l - start_ind == 0) {
                memcpy(&tmp, tmpin, 4);
            }
            else {
                tmpin += 4;
                memcpy(&tmp, tmpin, 4);
            }
            return tmp;
        }
        memcpy(&theta0, tmpin, sizeof(theta0));
        tmpin += sizeof(theta0);
        memcpy(&theta1, tmpin, sizeof(theta1));
        tmpin += sizeof(theta1);
        //std::cout<< "indexing "<<l<<std::endl;

        if (maxerror) {
            if (maxerror == 32) {
                tmp = read_bit_default(tmpin, maxerror, l - start_ind, theta1, theta0, maxerror);
            }
            else {
                tmp = read_bit_fix_T(tmpin ,maxerror, l - start_ind, theta1, theta0, 0);

            }
        }
        else {
            tmp = (theta0 + theta1 * (double)(l - start_ind));
        }

        // tmp = read_bit(tmpin ,maxerror , l-start_ind,theta1,theta0,0);
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
        std::cout << "Total block num is " << block_start_vec.size() << std::endl;
        return total_byte;
    }

    void destroy() {
        for (int i = 0;i < (int)block_start_vec.size();i++) {
            // block_start_vec[i].reset();
            free(block_start_vec[i]);
        }

    }
    std::string name() const {
        return "piecewise_fiting";
    }

};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
