
#ifndef PIECEWISE_COST_INTEGER_MERGE_HC_TEMPLATE_H_
#define PIECEWISE_COST_INTEGER_MERGE_HC_TEMPLATE_H_

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
    class Leco_cost_merge_hc
    {
    public:

        std::vector<uint8_t*> block_start_vec;
        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        std::vector<uint8_t*> block_start_vec_total;
        std::vector<uint32_t> segment_index_total;
        std::vector<uint32_t> segment_length_total;


        uint64_t total_byte_total = 0;
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

        uint32_t lower_bound(uint64_t v, uint32_t len, std::vector<uint32_t>& index)
        {
            uint32_t m;
            uint32_t x = 0;
            uint32_t y = len - 1;
            while (x <= y)
            {

                m = x + (y - x) / 2;
                if (v < index[m]) y = m - 1;
                else x = m + 1;
            }
            return y;

        }

        void newsegment(uint32_t origin_index, uint32_t end_index) {

            if(origin_index == end_index){
                return newsegment_1(origin_index,origin_index);
            }
            if(end_index == origin_index+1){
                return newsegment_2(origin_index,end_index);
            }
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
            float final_slope = mylr.theta1;
            double theta0 = mylr.theta0;

            T final_max_error = 0;
            std::vector<bool> signvec;
            std::vector<T> delta_final;

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                T tmp_val;
                int128_t pred = round(theta0 + final_slope * (double)(j - origin_index));
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
                delta_final_max_bit = bits_int_T<T>(final_max_error) + 1;
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
            // std::cout<<delta_final_max_bit<<std::endl;


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

        uint8_t* encodeArray8_int(T* in, const size_t length, uint8_t* res, size_t nvalue) {
            array = in;
            int64_t* delta_first_layer = new int64_t[length];
            int64_t* delta_second_layer = new int64_t[length];
            for (int i = nvalue*block_size;i < nvalue*block_size+length - 1;i++) {
                int64_t delta_first = (int64_t)in[i+1] -  (int64_t)in[i];
                delta_first_layer[i-nvalue*block_size] = delta_first; 
            }
            

            std::vector<uint32_t> segment_index_new;
            for (int i = 0;i < length - 2;i++) {
                delta_second_layer[i] = delta_first_layer[i + 1] - delta_first_layer[i];
            }

            // std::cout << "********" << 0 << "********" << std::endl;
            // for (auto item : segment_index_new) {
            //     std::cout << item << " ";
            // }
            // std::cout << std::endl;
            double model_cost = log((double)overhead)+log(8.0);
            int64_t max_second_layer_delta = INT64_MIN;
            int start_index = nvalue*block_size;
            int end_index = start_index;
            segment_index_new.push_back(start_index);
            // int tolerant_length = 0;
            while(end_index < nvalue*block_size+length-1){
                if(abs(delta_second_layer[end_index - nvalue*block_size]) > max_second_layer_delta){
                    max_second_layer_delta = abs(delta_second_layer[end_index-nvalue*block_size]);
                    // tolerant_length = exp((model_cost - log((double)max_second_layer_delta))/2.);
                }
                int tolerant_length = exp((model_cost - log((double)max_second_layer_delta))/2.);
                if(end_index-start_index<tolerant_length ){
                    end_index++;
                }
                else{
                    start_index = end_index+2;
                    end_index = start_index;
                    segment_index_new.push_back(start_index);
                    max_second_layer_delta = INT64_MIN;

                }
            }
            

            // for(auto item:segment_index_new){
            //     std::cout<<item<<" ";
            // }
            total_byte = 0;
            int segment_total = segment_index_new.size();
            segment_index_new.push_back(nvalue * block_size+length);
            for (int i=0;i<segment_total;i++) {
                newsegment(segment_index_new[i], segment_index_new[i+1]-1);
            }
            
            int iter = 0;
            uint64_t cost_decline = total_byte;
            while(cost_decline>0){
                
                iter++;
                cost_decline = total_byte;
                merge(nvalue);
                // merge_both_direction(nvalue);
                
                double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
                std::cout << "try "<<iter<<" segment number "<<(int)block_start_vec.size()<<" resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (sizeof(T) * block_size * 1.0);
                if(cost_decline_percent<0.01){
                    break;
                }
                
            }
            
            double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
            std::cout <<"segment number "<<(int)block_start_vec.size()<<" resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
            
            for(auto item: block_start_vec){
                block_start_vec_total.push_back(item);
            }
            for(auto item: segment_index){
                segment_index_total.push_back(item);
            }
            for(auto item: segment_length){
                segment_length_total.push_back(item);
            }
            total_byte_total+=total_byte;
            block_start_vec.clear();
            segment_index.clear();
            segment_length.clear();
            total_byte = 0;






            return res;

        }

        void merge(int nvalue) {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = segment_index[0]; // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int newsegment_num = 0;
            int total_segments = block_start_vec.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;
            segment_index.push_back((nvalue+1)*block_size);
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
                    start_index = block_size*(nvalue+1);
                    segment_num++;
                    break;
                }
                uint32_t init_cost = segment_length[segment_num] + segment_length[segment_num + 1];
                uint32_t merge_cost = 0;
                newsegment(start_index, segment_index[segment_num + 2] - 1);
                merge_cost = segment_length[total_segments + newsegment_num];
                if (init_cost > merge_cost) { // merge the two segments
                    // if(start_index==199999979){
                    //     std::cout<<"hi"<<std::endl;
                    // }
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
                    // if(start_index==199999979){
                    //     std::cout<<"hi"<<std::endl;
                    // }
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

        void merge_both_direction(int nvalue) {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = segment_index[0];  // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int total_segments = block_start_vec.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;
            segment_index.push_back(block_size*(nvalue+1));
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
                    start_index = block_size*(nvalue+1);
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
            new_segment_index.push_back(block_size*(nvalue+1));
            for (int i = 0;i < segment_number;i++) {
                newsegment(new_segment_index[i], new_segment_index[i + 1] - 1);
                // std::cout<<i<<" / "<<segment_index.size()<<" "<<total_byte<<std::endl;
            }
            new_segment_index.pop_back();

            block_start_vec.swap(new_block_start_vec);
            segment_index.swap(new_segment_index);
            segment_length.swap(new_segment_length);
            // std::cout << total_byte << std::endl;

        }

        uint32_t* decodeArray8(uint8_t* in, const size_t length, uint32_t* out, size_t nvalue) {
            //start_index + bit + theta0 + theta1 + numbers + delta

            return out;
        }


        T randomdecodeArray8(uint8_t *in, int to_find, uint32_t *out, size_t nvalue){

            uint32_t length = segment_index_total.size();
            uint8_t* this_block = block_start_vec_total[lower_bound(to_find, length, segment_index_total)];

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

            float theta1;
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);
            

            if(maxerror){
                // tmp_val = read_bit_fix_int_float<T>(tmpin, maxerror, to_find-start_ind, theta1, theta0);
                tmp_val = read_bit_fix_int<T>(tmpin, maxerror, to_find - start_ind, (double)theta1, theta0);
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
            std::cout << "Total block num is " << block_start_vec_total.size() << std::endl;
            return total_byte_total;
        }


    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
