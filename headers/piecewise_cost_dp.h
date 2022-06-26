
#ifndef PIECEWISE_COST_DP_H_
#define PIECEWISE_COST_DP_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {

    class piecewiseDp : public IntegerCODEC {
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
        int total_seg = 0;
        

        uint64_t total_byte = 0;
        // int overhead = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint64_t)*4;//start_index + start_key + slope
        // int overhead = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint8_t);
        int overhead = 10;
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

        uint32_t newsegment(uint32_t origin_index, uint32_t end_index) {
            int length = end_index - origin_index + 1;
            if(length == 1){
                uint32_t seg_len = newsegment_1(origin_index, origin_index);
                return seg_len;
            }
            if(length == 2){
                uint32_t seg_len = newsegment_2(origin_index, end_index);
                return seg_len;
            }
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(uint64_t)*2);
            uint8_t* out = descriptor;
            std::vector<double> indexes;
            std::vector<double> keys;
            for (int j = origin_index;j <= end_index;j++) {
                indexes.emplace_back(j - origin_index);
                keys.emplace_back(array[j]);
            }

            lr mylr;
            mylr.caltheta(indexes, keys, length);
            float final_slope = mylr.theta1;
            float theta0 = mylr.theta0;


            uint32_t final_max_error = 0;
            int* delta_final = new int[end_index - origin_index + 1];

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                long long  pred = (long long)theta0 + (long long)(final_slope * (float)(j - origin_index));
                uint32_t tmp_error = abs(pred - array[j]);
                delta_final[j - origin_index] = array[j] - pred;
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
            return segment_size;
        }

        uint32_t newsegment_2(uint32_t origin_index, uint32_t end_index) {
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
            return segment_size;
        }

        uint32_t newsegment_1(uint32_t origin_index, uint32_t end_index) {
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
            return segment_size;
        }


        uint8_t* encodeArray8(uint32_t* in, const size_t length, uint8_t* res, size_t nvalue) {
            array = in;
            std::vector<std::vector<uint32_t> > segment_start_index(nvalue,std::vector<uint32_t>());
            std::vector<std::vector<uint32_t> > segment_cost(nvalue,std::vector<uint32_t>());
            for(int i=0;i<nvalue;i++){
                for(int j=0;j<=i;j++){ //[j,i]
                    segment_cost[j].emplace_back( newsegment(j,i));
                    // std::cout<<"("<<j<<" , "<<i<<") "<<segment_cost[j][i -j]<<std::endl;
                }
                // std::cout<<"segment "<<i<<std::endl;
            }
            std::vector<uint64_t> segment_cost_sum(nvalue);
            for(int i=0;i<nvalue;i++){
                segment_cost_sum[i] = 0;
            }
            for(int i=0;i<nvalue;i++){

                uint64_t tmp = (1<<63) - 1;
                uint32_t start_ind = 0;
                uint32_t tmp_cost = segment_cost[0][i];
                if(tmp_cost<tmp){
                    tmp = tmp_cost;
                    start_ind = 0;
                }
                for(int j=1;j<=i;j++){
                    tmp_cost = segment_cost[j][i-j]+segment_cost_sum[j-1];
                    if(tmp_cost<tmp){
                        tmp = tmp_cost;
                        start_ind = j-1;
                    }
                }
                // std::cout<<"("<<start_ind<<" , "<<i<<") "<<tmp<<std::endl;
                segment_cost_sum[i] = tmp;
                if(start_ind!=0){
                    for(auto item : segment_start_index[start_ind]){
                        segment_start_index[i].push_back(item);

                    }
                    segment_start_index[i].emplace_back(start_ind+1);
                    
                }
                else{
                    segment_start_index[i].push_back(0);
                }
                
                // for(auto item : segment_start_index[i]){
                //     std::cout<<item<<" ";
                // }
                // std::cout<<std::endl;

            }
            destroy();
            segment_length.clear();
            block_start_vec.clear();
            segment_index.clear();
            total_byte = 0;
            segment_start_index[nvalue-1].push_back(nvalue);
            for(int i=0;i<segment_start_index[nvalue-1].size()-1;i++){
                newsegment(segment_start_index[nvalue-1][i], segment_start_index[nvalue-1][i+1]-1);
            }
            uint64_t min_cost = total_byte;
            total_seg += segment_start_index[nvalue-1].size() - 1 ;
            double compressrate = (min_cost) * 100.0 / (4 * block_size * 1.0);
            std::cout<<"resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
            std::cout<<"totalsegment: "<<total_seg<<std::endl;

            for (auto item : segment_start_index){
                item.clear();
            }
            for (auto item : segment_cost){
                item.clear();
            }
       
            segment_start_index.clear();
            segment_cost.clear();
            segment_cost_sum.clear();
            

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
            float theta0;
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
            if( maxerror==32){
                tmp = read_bit_default(tmpin,maxerror, l - start_ind, theta1, theta0, maxerror);
            } else{
                tmp = read_bit_fix_float_T(tmpin, maxerror, l - start_ind, theta1, theta0, 0);
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
            // std::cout << "Total block num is " << block_start_vec.size() << std::endl;
            return total_byte;
        }


        void destroy() {
            for (int i = 0;i < (int)block_start_vec.size();i++) {
                // block_start_vec[i].reset();
                free(block_start_vec[i]);
            }

        }
        std::string name() const {
            return "piecewise_cost_dp";
        }

    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */