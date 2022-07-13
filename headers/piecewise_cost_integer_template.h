
#ifndef PIECEWISE_COST_INTEGER_TEMPLATE_H_
#define PIECEWISE_COST_INTEGER_TEMPLATE_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {
    template <typename T>
    class Leco_cost
    {
    public:
        
        std::vector<uint8_t*> block_start_vec;
        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        uint64_t total_byte = 0;
        int overhead = 0;
        T* array;
        int block_num;
        int block_size;

        //start_index + bit + theta0 + theta1 + numbers + delta
        void init(int blocks, int blocksize, uint64_t delta) {
            block_num = blocks;
            block_size = blocksize;
            overhead = delta; // add some punishing item

        }

        uint32_t lower_bound(uint64_t v, uint32_t len)
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
            // if(origin_index == 199999979){
            //     std::cout<<"merging 199999979"<<std::endl;
            // }
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(T)*4);
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

            T final_max_error = 0;
            std::vector<bool> signvec;
            std::vector<T> delta_final;

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                T tmp_val;
                int128_t pred = round(theta0 + final_slope * (float)(j - origin_index));
                if ( array[j] > pred)
                {
                    tmp_val = array[j] - pred;
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = pred - array[j];
                    signvec.emplace_back(false); // means negative
                }

                delta_final.emplace_back(tmp_val);

                if (tmp_val > final_max_error)
                {
                    final_max_error = tmp_val;
                }
            }

            uint32_t delta_final_max_bit = 0;
            if(final_max_error){
                delta_final_max_bit=bits_int_T<T>(final_max_error) + 1;
            }
            
            
            if (delta_final_max_bit>= sizeof(T)*8){
                delta_final_max_bit = sizeof(T)*8;
            }

            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)delta_final_max_bit;
            out++;

            if(delta_final_max_bit== sizeof(T)*8){
                for (auto i = origin_index; i<=end_index; i++)
                {
                    memcpy(out, &array[i], sizeof(T));
                    out += sizeof(T);
                }
                uint64_t segment_size = out - descriptor;
                descriptor = (uint8_t*)realloc(descriptor, segment_size);
                block_start_vec.push_back(descriptor);
                segment_index.push_back(origin_index);
                segment_length.push_back(segment_size);
                total_byte += segment_size;
                
                return;
            }

            memcpy(out, &theta0, sizeof(theta0));
            out += sizeof(theta0);

            memcpy(out, &final_slope, sizeof(final_slope));
            out += sizeof(final_slope);


            if(delta_final_max_bit){
                out = write_delta_int_T(delta_final,signvec, out, delta_final_max_bit, (end_index - origin_index + 1));
            }


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
            
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(T)*2);
            uint8_t* out = descriptor;
            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)126; // this means that this segment only has two points
            out++;
            memcpy(out, &array[origin_index], sizeof(T));
            out += sizeof(T);
            memcpy(out, &(array[origin_index + 1]), sizeof(T));
            out += sizeof(T);

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);
            segment_length.push_back(segment_size);

            total_byte += segment_size;
        }

        void newsegment_1(uint32_t origin_index, uint32_t end_index) {

            uint8_t* descriptor = (uint8_t*)malloc(10 * sizeof(T));
            uint8_t* out = descriptor;
            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)127; // this means that this segment only has one point
            out++;
            memcpy(out, &array[origin_index], sizeof(T));
            out += sizeof(T);

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_length.push_back(segment_size);
            segment_index.push_back(origin_index);

            total_byte += segment_size;
        }

        uint8_t *encodeArray8_int(T *in, const size_t length, uint8_t *res, size_t nvalue){
            array = in;
            std::vector<uint32_t> indexes;
            for (uint32_t i = 0; i < block_size; i++) {
                indexes.push_back(i);
            }
            float high_slope = (float)INF;
            float low_slope = 0.;
            int128_t origin_key = in[0];
            uint32_t origin_index = indexes[0];
            uint32_t end_index = indexes[0];
            int tmp_delta_bit = 0;
            T tmp_max_delta = 0;
            for (int i = 1; i < (long long)block_size; i++) {

                int128_t key = in[i];
                int id = indexes[i];
                float tmp_point_slope = ((key - origin_key) + 0.0) / ((id - origin_index) + 0.0);
                if (id == block_size - 1) {
                    if(id==origin_index){
                        newsegment_1(origin_index, origin_index);
                        break;
                    }
                    newsegment(origin_index, id);
                    break;
                }
                if(id==origin_index){
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
                        T pred = origin_key + (float)(id - origin_index) * tmp_slope;
                        T tmp_error = 0;
                        if(pred > in[j]){
                            tmp_error = pred - in[j];
                        }
                        else{
                            tmp_error = in[j] - pred;
                        }
                        if (tmp_error > tmp_max_delta) {
                            tmp_max_delta = tmp_error;
                        }
                        tmp_delta_bit = bits_int_T<T>(tmp_max_delta) + 1;
                    } 
                    int new_cost = tmp_delta_bit * (id - origin_index + 1);
                    int old_cost =  sizeof(double)+sizeof(double)+sizeof(uint32_t);
                    // if (tmp_max_delta >= tolerance) {
                    //     newsegment_2(origin_index, origin_index+1);

                    //     high_slope = (float)INF;
                    //     low_slope = 0.0;
                    //     origin_index = id;
                    //     origin_key = key;
                    //     end_index = id;
                    //     tmp_delta_bit = 0;
                    //     tmp_max_delta = 0;
                    //     continue;

                    // }
                    if(new_cost-old_cost>overhead){
                        newsegment_2(origin_index, origin_index+1);
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
                T pred = origin_key + (float)(id - origin_index) * tmp_slope;
                T tmp_error = 0;
                if(pred>key){
                    tmp_error = pred - key;
                }
                else{
                    tmp_error = key - pred;
                }
                
                int tmp_error_bit = bits_int_T<T>(tmp_error) + 1;
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
                    T pred = origin_key + (float)(id - origin_index) * mid_slope;
                    T tmp_error = 0;
                    if(pred>in[id]){
                        tmp_error = pred - in[id];
                    }
                    else{
                        tmp_error = in[id] - pred;
                    }
                    int delta_max_bit = bits_int_T<T>(tmp_error) + 1;
                    uint64_t cost = (id - origin_index + 1) * (delta_max_bit - tmp_delta_bit);
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
                        if (id == block_size - 1) {
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
            while(cost_decline>0){
                
                iter++;
                cost_decline = total_byte;
                merge();
                
                double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
                std::cout << "try "<<iter<<" segment number "<<(int)block_start_vec.size()<<" resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (sizeof(T) * block_size * 1.0);
                if(cost_decline_percent<0.01){
                    break;
                }
                
            }
            // for(int j=0;j<4;j++){
            //     merge();
            // }


            return res;

        }

        void merge(){
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
            while(segment_num < total_segments){
                // std::cout<<"segment_num: "<<segment_num <<" / "<<total_segments<<std::endl;

                if (segment_num == total_segments - 1) {
                    // std::cout <<segment_num<<"///"<<total_segments<<" "<< block_start_vec[segment_num] << std::endl;
                    new_block_start_vec.push_back(block_start_vec[segment_num]);
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index=block_size;
                    segment_num++;
                    break;
                }
                uint32_t init_cost = segment_length[segment_num] + segment_length[segment_num+1];
                uint32_t merge_cost = 0;
                newsegment(start_index, segment_index[segment_num+2]-1);
                merge_cost = segment_length[total_segments + newsegment_num];
                if(init_cost>merge_cost){ // merge the two segments
                    // if(start_index==199999979){
                    //     std::cout<<"hi"<<std::endl;
                    // }
                    new_block_start_vec.emplace_back(block_start_vec[total_segments+newsegment_num]);
                    new_segment_index.emplace_back(start_index);
                    new_segment_length.emplace_back(merge_cost);
                    totalbyte_after_merge += merge_cost;
                    start_index=segment_index[segment_num+2];
                    
                    segment_num+=2;
                    // std::cout<<segment_num<<std::endl;
                    newsegment_num++;
                }
                else {
                    // if(start_index==199999979){
                    //     std::cout<<"hi"<<std::endl;
                    // }
                    // std::cout <<segment_num<<"/"<<total_segments<<" "<< block_start_vec[segment_num] << std::endl;
                    new_block_start_vec.emplace_back(block_start_vec[segment_num]);
                    // new_block_start_vec.emplace_back(std::move(std::unique_ptr<uint8_t>(block_start_vec[segment_num])));
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index=segment_index[segment_num+1];
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

        uint32_t* decodeArray8(uint8_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            //start_index + bit + theta0 + theta1 + numbers + delta

            return out;
        }


        T randomdecodeArray8(uint8_t *in, int to_find, uint32_t *out, size_t nvalue){

            uint32_t length = segment_index.size();
            uint8_t* this_block = block_start_vec[lower_bound(to_find, length)];

            uint8_t* tmpin = this_block;
            
            uint32_t start_ind;
            memcpy(&start_ind, tmpin, 4);
            tmpin += 4;

            uint8_t maxerror;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            if(maxerror==sizeof(T)*8){
                T tmp_val = reinterpret_cast<T *>(tmpin)[to_find];
                return tmp_val;
            }

            T tmp_val = 0;
            if (maxerror == 127) {
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                return tmp_val;
            }
            if (maxerror == 126) {
                if (to_find - start_ind == 0) {
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                }
                else {
                    tmpin += sizeof(tmp_val);
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                }
                return tmp_val;
            }

            double theta0;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);

            double theta1;
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);
            

            if(maxerror){
                // tmp_val = read_bit_fix_int_float<T>(tmpin, maxerror, to_find-start_ind, theta1, theta0);
                tmp_val = read_bit_fix_int<T>(tmpin, maxerror, to_find - start_ind, theta1, theta0);
            } 
            else{
                tmp_val = (T)round(theta0+theta1 * (double)(to_find - start_ind));
            }
  
            
                
            return tmp_val;

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


    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
