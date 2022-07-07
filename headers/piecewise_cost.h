
#ifndef PIECEWISE_COST_H_
#define PIECEWISE_COST_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {

    class piecewiseCost : public IntegerCODEC {
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
        // int overhead = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint64_t)*4;//start_index + start_key + slope
        // int overhead = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint8_t);
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
            float final_slope = mylr.theta1;
            double theta0 = mylr.theta0;

            uint32_t final_max_error = 0;
            int* delta_final = new int[end_index - origin_index + 1];

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                long long  pred = round(theta0 + final_slope * (double)(j - origin_index));

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


            // if(1636>=origin_index && 1636<=end_index){
                // std::cout<<"("<<origin_index<<" , "<<end_index<<") "<<(end_index - origin_index + 1)<<" "<<array[origin_index]<<" "<<array[end_index]<<" "<<delta_final_max_bit<<std::endl;
                // for(int j=origin_index;j<=end_index;j++){
                //     std::cout<<delta_final[j - origin_index]<<" ";
                // }
                // std::cout<<std::endl;
            // }


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
                if (id == nvalue - 1) {
                    newsegment(origin_index, id);
                    break;
                }
                if (id == origin_index) {
                    continue;
                }
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
                        tmp_delta_bit = bits(tmp_max_delta) + 1;
                    }
                    int new_cost = tmp_delta_bit * (id - origin_index + 1);
                    int old_cost = sizeof(double) + sizeof(double) + sizeof(uint32_t) + sizeof(uint8_t);

                    if (new_cost - old_cost > overhead) {
                        // std::cout<<"newseg "<<origin_index<<" "<<id<<" "<<new_cost<<std::endl;
                        newsegment_2(origin_index, origin_index + 1);
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
                if (tmp_error_bit <= tmp_delta_bit) {
                    end_index = id;
                    if (tmp_error > tmp_max_delta) {
                        tmp_max_delta = tmp_error;
                    }
                    continue;
                }
                else {
                    float mid_slope = (high_slope + low_slope) / 2.;
                    if (tmp_point_slope < low_slope) {
                        mid_slope = (high_slope + tmp_point_slope) / 2.;
                    }
                    if (tmp_point_slope > high_slope) {
                        mid_slope = (low_slope + tmp_point_slope) / 2.;
                    }
                    long long pred = origin_key + (float)(id - origin_index) * mid_slope;
                    int tmp_error = abs(pred - in[id]);
                    int delta_max_bit = bits(tmp_error) + 1;
                    int cost = (id - origin_index + 1) * (delta_max_bit - tmp_delta_bit);
                    // std::cout<<origin_index<<" "<<id<<" "<<cost<<std::endl;
                    if (cost < overhead) {

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


                    }
                    else {
                        // delete[] delta;
                        // write the last segment & start a new segment
                        newsegment(origin_index, end_index);
                        if (id == nvalue - 1) {
                            newsegment_1(id, id);
                        }
                        high_slope = (float)INF;
                        low_slope = 0.0;
                        origin_index = id;
                        origin_key = key;
                        end_index = id;

                        tmp_delta_bit = 0;
                        tmp_max_delta = 0;


                    }

                }

            }
            int iter = 0;
            uint64_t cost_decline = total_byte;
            while (cost_decline > 0) {

                iter++;
                cost_decline = total_byte;
                // merge();
                merge_both_direction();

                double compressrate = (total_byte) * 100.0 / (4 * block_size * 1.0);
                std::cout << "try " << iter << " segment number " << (int)block_start_vec.size() << " resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (4 * block_size * 1.0);
                if (cost_decline_percent < 0.01) {
                    break;
                }

            }
            // for (auto item : segment_index) {
            //     std::cout << item << std::endl;
            // }


            return res;

        }

        void merge() {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = 0; // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int newsegment_num = 0;
            int total_segments = block_start_vec.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;
            segment_index.push_back(block_size);
            std::vector<uint8_t*> new_block_start_vec;
            std::vector<uint32_t> new_segment_index;
            std::vector<uint32_t> new_segment_length;
            while (segment_num < total_segments) {
                // std::cout<<"segment_num: "<<segment_num <<" / "<<total_segments<<std::endl;

                if (segment_num == total_segments - 1) {
                    // std::cout <<segment_num<<"///"<<total_segments<<" "<< block_start_vec[segment_num] << std::endl;
                    new_block_start_vec.push_back(block_start_vec[segment_num]);
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = block_size;
                    segment_num++;
                    break;
                }
                uint32_t init_cost = segment_length[segment_num] + segment_length[segment_num + 1];
                uint32_t merge_cost = 0;
                newsegment(start_index, segment_index[segment_num + 2] - 1);
                merge_cost = segment_length[total_segments + newsegment_num];
                if (init_cost > merge_cost) { // merge the two segments
                    // new_block_start_vec.emplace_back(std::unique_ptr<uint8_t>(block_start_vec[total_segments+newsegment_num]));
                    // if(segment_num ==1243665 ){
                    //     std::cout <<segment_num<<"//"<<total_segments<<" "<< block_start_vec[total_segments+newsegment_num] << std::endl;
                    // }
                    // std::cout<<"merge "<<segment_num<<" "<<segment_num+1<<" ( "<<start_index<<" , "<<segment_index[segment_num+2]-1<<" ) "<<" init cost: "<<init_cost<<" merge cost: "<<merge_cost<<std::endl;

                    new_block_start_vec.emplace_back(block_start_vec[total_segments + newsegment_num]);
                    new_segment_index.emplace_back(start_index);
                    new_segment_length.emplace_back(merge_cost);
                    totalbyte_after_merge += merge_cost;
                    start_index = segment_index[segment_num + 2];
                    segment_num += 2;
                    // std::cout<<segment_num<<std::endl;
                    newsegment_num++;
                }
                else {
                    // std::cout <<segment_num<<"/"<<total_segments<<" "<< block_start_vec[segment_num] << std::endl;
                    new_block_start_vec.emplace_back(block_start_vec[segment_num]);
                    // new_block_start_vec.emplace_back(std::move(std::unique_ptr<uint8_t>(block_start_vec[segment_num])));
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = segment_index[segment_num + 1];
                    segment_num++;
                    newsegment_num++;
                }

            }
            block_start_vec.swap(new_block_start_vec);
            segment_index.swap(new_segment_index);
            segment_length.swap(new_segment_length);
            total_byte = totalbyte_after_merge;
            // std::cout<<total_byte<<std::endl;

        }


        void merge_both_direction() {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = 0; // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int total_segments = block_start_vec.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;
            segment_index.push_back(block_size);
            std::vector<uint8_t*> new_block_start_vec;
            std::vector<uint32_t> new_segment_index;
            std::vector<uint32_t> new_segment_length;
            new_block_start_vec.push_back(block_start_vec[segment_num]);
            new_segment_index.emplace_back(segment_index[segment_num]);
            new_segment_length.emplace_back(segment_length[segment_num]);
            totalbyte_after_merge += segment_length[segment_num];
            segment_num++;

            while (segment_num < total_segments) {
                // std::cout<<"segment_num: "<<segment_num <<" / "<<total_segments<<std::endl;

                if (segment_num == total_segments - 1) {
                    // only can try merging with former one
                    new_block_start_vec.push_back(block_start_vec[segment_num]);
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = block_size;
                    segment_num++;
                    break;
                }
                int last_merged_segment = new_segment_length.size() - 1;

                uint32_t init_cost_front = segment_length[segment_num] + new_segment_length[last_merged_segment];
                newsegment(new_segment_index[last_merged_segment], segment_index[segment_num + 1] - 1);
                uint32_t merge_cost_front = segment_length[segment_length.size() - 1];
                int saved_cost_front = init_cost_front - merge_cost_front;

                uint32_t init_cost_back = segment_length[segment_num] + segment_length[segment_num + 1];
                newsegment(segment_index[segment_num], segment_index[segment_num + 2] - 1);
                uint32_t merge_cost_back = segment_length[segment_length.size() - 1];
                int saved_cost_back = init_cost_back - merge_cost_back;

                int saved_cost = std::max(saved_cost_front, saved_cost_back);
                if (saved_cost <= 0) {
                    // do not merge
                    new_block_start_vec.emplace_back(block_start_vec[segment_num]);
                    // new_block_start_vec.emplace_back(std::move(std::unique_ptr<uint8_t>(block_start_vec[segment_num])));
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    // std::cout<<"not merge "<<totalbyte_after_merge<<std::endl;
                    start_index = segment_index[segment_num + 1];
                    segment_num++;

                    continue;
                }
                if (saved_cost_back > saved_cost_front) {
                    // merge with back
                    new_block_start_vec.emplace_back(block_start_vec[block_start_vec.size() - 1]);
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(merge_cost_back);
                    totalbyte_after_merge += merge_cost_back;
                    // std::cout<<"merge with back "<<totalbyte_after_merge<<std::endl;
                    start_index = segment_index[segment_num + 2];
                    segment_num += 2;
                    // std::cout<<segment_num<<std::endl;

                }
                else {
                    // merge with front
                    new_block_start_vec[new_block_start_vec.size() - 1] = block_start_vec[block_start_vec.size() - 2];
                    totalbyte_after_merge -= new_segment_length[new_segment_length.size() - 1];
                    new_segment_length[new_segment_length.size() - 1] = merge_cost_front;
                    totalbyte_after_merge += merge_cost_front;
                    // std::cout<<"merge with front "<<totalbyte_after_merge<<std::endl;
                    start_index = segment_index[segment_num + 1];
                    segment_num += 1;
                    // std::cout<<segment_num<<std::endl;

                }

            }
            total_byte = 0;
            block_start_vec.clear();
            segment_index.clear();
            segment_length.clear();
            int segment_number = (int)new_segment_index.size();
            new_segment_index.push_back(block_size);
            for (int i = 0;i < segment_number;i++) {
                newsegment(new_segment_index[i], new_segment_index[i + 1] - 1);
                // std::cout<<i<<" / "<<segment_index.size()<<" "<<total_byte<<std::endl;
            }
            new_segment_index.pop_back();

            block_start_vec.swap(new_block_start_vec);
            segment_index.swap(new_segment_index);
            segment_length.swap(new_segment_length);
            std::cout << total_byte << std::endl;

        }

        uint32_t* decodeArray8(uint8_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            //start_index + bit + theta0 + theta1 + numbers + delta

            return out;
        }


        uint32_t randomdecodeArray8(uint8_t* in, const size_t l, uint32_t* out, size_t nvalue) {

            uint32_t length = segment_index.size();
            uint8_t* this_block = block_start_vec[lower_bound(l, length)];

            uint8_t* tmpin = this_block;
            double theta0;
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
                    tmp = read_bit_fix_int<uint32_t>(tmpin, maxerror, l - start_ind, (double)theta1, theta0);

                }
            }
            else {
                tmp = round(theta0 + theta1 * (double)(l - start_ind));
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
            return "piecewise_cost";
        }

    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
