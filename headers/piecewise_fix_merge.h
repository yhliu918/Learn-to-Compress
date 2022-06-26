
#ifndef PIECEWISE_FIX_MERGE_H_
#define PIECEWISE_FIX_MERGE_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {

    class piecewise_fix_merge : public IntegerCODEC {
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
        // int overhead = 18;
        uint32_t* array;
        int tolerance = 0;
        int block_num;
        int block_size;

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
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t)*2);
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
                long long  pred =  (long long)(theta0 + final_slope * (double)(j - origin_index));
                uint32_t tmp_error = abs(pred - (long long)array[j]);
                delta_final[j - origin_index] = (long long)array[j] - pred;
                if (tmp_error > final_max_error) {
                    final_max_error = tmp_error;
                }
            }
            uint32_t delta_final_max_bit = bits_int_T<uint32_t>(final_max_error) + 1;

            
            
            if (delta_final_max_bit>=32){
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
            if(delta_final_max_bit==32){
                out = write_delta_default(array+origin_index, out,delta_final_max_bit, end_index - origin_index + 1);
            }
            else{
                out = write_delta_T(delta_final, out, delta_final_max_bit, (end_index - origin_index + 1));
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


        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* res, size_t nvalue) {
            // nvalue is block_number here
            array = in;
            int block_interval = block_size / nvalue;
            nvalue = block_size / block_interval;
            if(block_size > block_interval * nvalue){
                nvalue++;
            }
            for(int i=0;i<nvalue;i++){
                int start_index = block_interval*i;
                int end_index = std::min(block_interval*(i+1), block_size) - 1;
                // std::cout<<"block "<<i<<" "<<start_index<<" "<<end_index<<std::endl;
                newsegment(start_index, end_index);
            }

            int iter = 0;
            double compressrate = (total_byte) * 100.0 / (4 * block_size * 1.0);
            std::cout << "try "<<iter<<" segment number "<<(int)block_start_vec.size()<<" resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;

            uint64_t cost_decline = total_byte;
            while(cost_decline>0){
                
                iter++;
                cost_decline = total_byte;
                merge();
                
                compressrate = (total_byte) * 100.0 / (4 * block_size * 1.0);
                std::cout << "try "<<iter<<" segment number "<<(int)block_start_vec.size()<<" resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (4 * block_size * 1.0);
                if(cost_decline_percent<0.01){
                    break;
                }
                
            }
            
            // merge();

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
                    // new_block_start_vec.emplace_back(std::unique_ptr<uint8_t>(block_start_vec[total_segments+newsegment_num]));
                    // std::cout<<"merge "<<segment_num<<" "<<segment_num+1<<" ( "<<start_index<<" , "<<segment_index[segment_num+2]-1<<" ) "<<" init cost: "<<init_cost<<" merge cost: "<<merge_cost<<std::endl;
                    
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
            if( maxerror==32){
                tmp = read_bit_default(tmpin,maxerror, l - start_ind, theta1, theta0, maxerror);
            } else{
                tmp = read_bit_fix_T(tmpin, maxerror, l - start_ind, theta1, theta0, 0);
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
            return "piecewise_fix_merge";
        }

    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
