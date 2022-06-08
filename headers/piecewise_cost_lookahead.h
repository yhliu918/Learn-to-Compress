
#ifndef PIECEWISE_COST_AHEAD_H_
#define PIECEWISE_COST_AHEAD_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#define INF 0x7f7fffff

namespace Codecset {

    class piecewiseCostAhead : public IntegerCODEC {
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
        uint32_t counter = 0;
        uint32_t total_byte = 0;
        int overhead = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint64_t) * 4;//start_index + start_key + slope
        // int overhead = 36;//start_index + start_key + slope
        int LookAheadLen = 20;
        int tolerance = 0;
        int block_num;
        int block_size;
        uint32_t* array;

        //start_index + bit + theta0 + theta1 + numbers + delta
        void init(int blocks, int blocksize, int delta) {
            block_num = blocks;
            block_size = blocksize;
            tolerance = delta; // add some punishing item

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
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t));
            uint8_t* out = descriptor;
            int length = end_index - origin_index + 1;
            double* indexes = new double[length];
            double* keys = new double[length];
            for (int j = origin_index;j <= end_index;j++) {
                indexes[j - origin_index] = j - origin_index;
                keys[j - origin_index] = array[j];
            }

            lr mylr;
            mylr.caltheta(indexes, keys, length);
            float final_slope = mylr.theta1;
            int32_t theta0 = mylr.theta0;

            int final_max_error = 0;
            int* delta_final = new int[end_index - origin_index + 1];

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                long long  pred = (long long)theta0 + (long long)(final_slope * (float)(j - origin_index));
                int tmp_error = abs(pred - array[j]);
                delta_final[j - origin_index] = array[j] - pred;
                if (tmp_error > final_max_error) {
                    final_max_error = tmp_error;
                }
            }
            int delta_final_max_bit = bits(final_max_error) + 1;
            // if (end_index<2100){
            //     std::cout<<origin_index<<" "<<end_index<<" "<<(long long)theta0<<" "<<final_slope<<std::endl;
            // }


            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)delta_final_max_bit;
            out++;

            memcpy(out, &theta0, sizeof(theta0));
            out += sizeof(theta0);

            memcpy(out, &final_slope, sizeof(final_slope));
            out += sizeof(final_slope);
            out = write_delta_T(delta_final, out, delta_final_max_bit, (end_index - origin_index + 1));


            // if(3530>=origin_index && 3530<=end_index){
            //     std::cout<<"("<<origin_index<<" , "<<end_index<<") "<<(end_index - origin_index + 1)<<" "<<in[origin_index]<<" "<<in[end_index]<<" "<<delta_final_max_bit<<" "<<(long long)theta0<<" "<<final_slope<<std::endl;
            //     for(int j=origin_index;j<=end_index;j++){
            //         std::cout<<delta_final[j - origin_index]<<" ";
            //     }
            //     std::cout<<std::endl;
            // }


            delete[] delta_final;


            int segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);
            total_byte += segment_size;
            // if(origin_index == 2024){
            //     std::cout<<segment_size<<" "<<end_index<<std::endl;
            // }
        }

        void newsegment_2(uint32_t origin_index, uint32_t end_index) {
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

            int segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);

            total_byte += segment_size;
        }

        void newsegment_1(uint32_t origin_index, uint32_t end_index) {
            uint8_t* descriptor = (uint8_t*)malloc(3 * sizeof(uint64_t));
            uint8_t* out = descriptor;
            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)127; // this means that this segment only has one point
            out++;
            memcpy(out, &array[origin_index], sizeof(uint32_t));
            out += sizeof(uint32_t);

            int segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);

            total_byte += segment_size;
        }

        uint32_t lookahead(uint32_t origin_index, uint32_t lookaheadlength) {
            // return total byte of [origin_index, origin_index + lookaheadlength]
            uint32_t end_index = origin_index + lookaheadlength;
            if (block_size-1 < end_index) {
                end_index = block_size - 1;
            }
            float high_slope = (float)INF;
            float low_slope = 0.;
            if (origin_index + 1 <= end_index) {
                float tmp_point_slope = ((array[origin_index + 1] - array[origin_index]) + 0.0) / (1. + 0.0);
                low_slope = tmp_point_slope;
            }
            else {
                return 9;
            }
            if (origin_index + 2 <= end_index) {
                float tmp = 0;
                float tmp_point_slope = ((array[origin_index + 2] - array[origin_index]) + 0.0) / (2. + 0.0);
                if (tmp_point_slope < low_slope) {
                    tmp = low_slope;
                    low_slope = tmp_point_slope;
                    high_slope = tmp;
                }
                else {
                    high_slope = tmp_point_slope;
                }
            }
            else {
                return 13;
            }

            uint32_t totalcost = 13;
            for (uint32_t i = origin_index + 3;i <= end_index;i++) {
                float tmp_point_slope = ((array[i] - array[origin_index]) + 0.0) / (i - origin_index + 0.0);
                if (tmp_point_slope < low_slope) {
                    low_slope = tmp_point_slope;
                }
                else if (tmp_point_slope > high_slope) {
                    high_slope = tmp_point_slope;
                }
            }
            float mid_slope = (low_slope + high_slope) / 2.;
            uint32_t final_max_error = 0;
            for (uint32_t i = origin_index + 1;i <= end_index;i++) {
                long long  pred = (long long)array[origin_index] + (long long)(mid_slope * (float)(i - origin_index));
                int tmp_error = abs(pred - array[i]);
                if (tmp_error > final_max_error) {
                    final_max_error = tmp_error;
                }
            }
            int delta_final_max_bit = bits(final_max_error) + 1;
            totalcost += ceil(delta_final_max_bit * lookaheadlength/8.);
            return totalcost;


        }

        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* res, size_t nvalue) {
            array = in;
            std::vector<uint32_t> indexes;
            for (uint32_t i = 0; i < nvalue; i++) {
                indexes.push_back(i);
            }
            float high_slope = (float)INF;
            float low_slope = 0.;
            long long origin_key = in[0];
            int origin_index = indexes[0];
            int end_index = indexes[0];
            int tmp_delta_bit = 0;
            int tmp_max_delta = 0;
            for (int i = 1; i < (long long)nvalue; i++) {
                long long key = in[i];
                int id = indexes[i];
                float tmp_point_slope = ((key - origin_key) + 0.0) / ((id - origin_index) + 0.0);
                if (id == origin_index + 1) {
                    low_slope = tmp_point_slope;
                    end_index = id;
                    continue;
                }
                if (id == origin_index + 2) {
                    float tmp = 0;
                    if (tmp_point_slope < low_slope) {
                        tmp = low_slope;
                        low_slope = tmp_point_slope;
                        high_slope = tmp;
                    }
                    else {
                        high_slope = tmp_point_slope;
                    }
                    end_index = id;
                    float tmp_slope = (high_slope + low_slope) / 2;
                    for (int j = origin_index + 1;j < id;j++) {
                        long long pred = origin_key + (float)(id - origin_index) * tmp_slope;
                        int tmp_error = abs(pred - in[j]);
                        if (tmp_error > tmp_max_delta) {
                            tmp_max_delta = tmp_error;
                        }
                    }
                    tmp_delta_bit = bits(tmp_max_delta) + 1;
                    double pre_ratio = total_byte / (double)(origin_index * sizeof(uint32_t));
                    double tmp_ratio = tmp_delta_bit / 32.0;
                    // std::cout<<"pre_ratio: "<<pre_ratio<<" tmp_ratio: "<<tmp_ratio<<std::endl;

                    if (total_byte != 0 && tmp_ratio > pre_ratio) {

                        newsegment_2(origin_index, end_index);
                        high_slope = (float)INF;
                        low_slope = 0.0;
                        origin_index = id;
                        origin_key = key;
                        end_index = id;
                        tmp_delta_bit = 0;
                        tmp_max_delta = 0;
                        continue;

                    }

                    continue;
                }

                float tmp_slope = (high_slope + low_slope) / 2;
                long long pred = origin_key + (float)(id - origin_index) * tmp_slope;
                int tmp_error = abs(pred - key);
                int tmp_error_bit = bits(tmp_error) + 1;
                // check whether this point is assigned to the tmp segment
                // by calculating the increase of the error this point brings to the segment
                if (tmp_error_bit <= tmp_delta_bit) {
                    if (id == nvalue - 1) {
                        newsegment(origin_index, id);
                        break;
                    }
                    end_index = id;
                    if (tmp_error > tmp_max_delta) {
                        tmp_max_delta = tmp_error;
                    }
                    continue;
                }
                else {
                    uint32_t totalbytecopy = total_byte;
                    float mid_slope = (high_slope + low_slope) / 2.;
                    if (tmp_point_slope < low_slope) {
                        mid_slope = (high_slope + tmp_point_slope) / 2.;
                    }
                    if (tmp_point_slope > high_slope) {
                        mid_slope = (low_slope + tmp_point_slope) / 2.;
                    }

                    // Choice A: with current segment
                    long long pred = origin_key + (float)(id - origin_index) * mid_slope;
                    int tmp_error = abs(pred - in[id]);
                    int delta_max_bit = bits(tmp_error) + 1;
                    int cost = 13 + ceil(((id - origin_index + 1) * std::max(tmp_delta_bit, delta_max_bit))/8.);
                    double currentsegmentCR = (double)(totalbytecopy + cost) / (double)((id+1) * sizeof(uint32_t));

                    // Choice B: with new segment
                    uint32_t lookahead_cost = lookahead(id, LookAheadLen);
                    int currentsegcost = 13 + ceil(((id - origin_index) * tmp_delta_bit)/8.);
                    double newsegmentCR = (double)(totalbytecopy + currentsegcost + lookahead_cost) / (double)((id + LookAheadLen+1) * sizeof(uint32_t));

                    // Choice C: outlier
                    int outlier_cost = 13;
                    double outlierCR = (double)(totalbytecopy + currentsegcost + outlier_cost) / (double)((id + 2) * sizeof(uint32_t));

                    // choose the one with the lowest CR
                    if (currentsegmentCR < newsegmentCR && currentsegmentCR < outlierCR) {
                        if (tmp_point_slope < low_slope) {
                            low_slope = tmp_point_slope;
                        }
                        if (low_slope < 0) {
                            low_slope = 0.0;
                        }
                        if (tmp_point_slope > high_slope) {
                            high_slope = tmp_point_slope;
                        }
                        end_index = id;
                        if (delta_max_bit > tmp_delta_bit) {
                            tmp_delta_bit = delta_max_bit;
                        }
                        continue;
                    }
                    else if (newsegmentCR < outlierCR && newsegmentCR < currentsegmentCR) {
                        newsegment(origin_index, end_index);
                        if (id == nvalue - 1) {
                            newsegment_1(id, id);
                            high_slope = (float)INF;

                        }
                        high_slope = (float)INF;
                        low_slope = 0.0;
                        origin_index = id;
                        origin_key = key;
                        end_index = id;
                        tmp_delta_bit = 0;
                        tmp_max_delta = 0;
                        continue;


                    }
                    else if (outlierCR <= newsegmentCR && outlierCR < currentsegmentCR) {
                        newsegment(origin_index, end_index);
                        newsegment_2(id, id + 1);
                        high_slope = (float)INF;
                        low_slope = 0.0;
                        origin_index = id + 1;
                        origin_key = in[id + 1];
                        end_index = id + 1;
                        tmp_delta_bit = 0;
                        tmp_max_delta = 0;
                        continue;

                    }



                }

            }

            return res;

        }


        uint32_t* decodeArray8(uint8_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            //start_index + bit + theta0 + theta1 + numbers + delta

            return out;
        }


        uint32_t randomdecodeArray8(uint8_t* in, const size_t l, uint32_t* out, size_t nvalue) {

            uint32_t length = segment_index.size();
            uint8_t* this_block = block_start_vec[lower_bound(l, length)];

            uint8_t* tmpin = this_block;
            int32_t theta0;
            float theta1;
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
            memcpy(&theta0, tmpin, 4);
            tmpin += 4;
            memcpy(&theta1, tmpin, 4);
            tmpin += 4;
            //std::cout<< "indexing "<<l<<std::endl;
            tmp = read_bit_fix_float_T(tmpin, maxerror, l - start_ind, theta1, theta0, 0);
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
                free(block_start_vec[i]);
            }

        }
        std::string name() const {
            return "piecewise_cost_ahead";
        }

    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
