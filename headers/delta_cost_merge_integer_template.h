
#ifndef DELTA_COST_INTEGER_MERGE_TEMPLATE_H_
#define DELTA_COST_INTEGER_MERGE_TEMPLATE_H_

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
    class Delta_cost_merge
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
        uint32_t cal_length(uint32_t origin_index, uint32_t end_index) {
            uint32_t totallength = 0;
            T max_error = 0;
            for (auto i = origin_index; i <= end_index - 1; i++)
            {
                T tmp_val;
                if (array[i + 1] > array[i])
                {
                    tmp_val = array[i + 1] - array[i];
                }
                else
                {
                    tmp_val = array[i] - array[i + 1];
                }

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }

            uint8_t max_bit = 0;
            if (max_error)
            {
                max_bit = bits_int_T(max_error) + 1;
            }
            totallength += (sizeof(uint32_t) + sizeof(uint8_t) + sizeof(T));
            totallength += ceil(max_bit * (end_index - origin_index + 1) / 8);
            return totallength;

        }

        void newsegment(uint32_t origin_index, uint32_t end_index) {

            if (origin_index == end_index) {
                return newsegment_1(origin_index, origin_index);
            }
            if (origin_index == end_index + 1) {
                return newsegment_2(origin_index, origin_index);
            }
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(T) * 8);
            uint8_t* out = descriptor;

            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;

            for (auto i = origin_index; i <= end_index - 1; i++)
            {
                T tmp_val;
                if (array[i + 1] > array[i])
                {
                    tmp_val = array[i + 1] - array[i];
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = array[i] - array[i + 1];
                    signvec.emplace_back(false); // means negative
                }

                delta.emplace_back(tmp_val);

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }



            uint8_t max_bit = 0;
            if (max_error)
            {
                max_bit = bits_int_T(max_error) + 1;
            }

            if (max_bit > sizeof(T) * 8) {
                max_bit = sizeof(T) * 8;
            }
            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            memcpy(out, &max_bit, sizeof(max_bit));
            out += sizeof(max_bit);
            if (max_bit == sizeof(T) * 8) {
                for (auto i = origin_index; i <= end_index; i++)
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

            memcpy(out, &array[origin_index], sizeof(T));
            out += sizeof(T);

            if (max_bit)
            {
                out = write_delta_int_T(delta, signvec, out, max_bit, (end_index - origin_index));
                // because delta has N-1 deltas to write
            }


            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec.push_back(descriptor);
            segment_index.push_back(origin_index);
            segment_length.push_back(segment_size);
            total_byte += segment_size;
            // std::cout<<"segment_size: "<<segment_size<<std::endl;
            // if(origin_index == 2024){
            //     std::cout<<segment_size<<" "<<end_index<<std::endl;
            // }
        }

        void newsegment_2(uint32_t origin_index, uint32_t end_index) {

            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(T) * 2);
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
            int* delta_first_layer = new int[length];
            int* delta_second_layer = new int[length];
            for (int i = nvalue * block_size;i < nvalue * block_size + length - 1;i++) {
                T delta_first = 0;
                if (in[i + 1] > in[i]) {
                    delta_first = in[i + 1] - in[i];
                }
                else {
                    delta_first = in[i] - in[i + 1];
                }
                if (delta_first) {
                    delta_first_layer[i - nvalue * block_size] = bits_int_T<T>(delta_first) + 1;
                }
                else {
                    delta_first_layer[i - nvalue * block_size] = 0;
                }
            }
            int max_second_delta = 0;
            // max int of int128
            int min_second_delta = (1 << 30);
            for (int i = 0;i < length - 2;i++) {
                delta_second_layer[i] = delta_first_layer[i + 1] - delta_first_layer[i];
                if (abs(delta_second_layer[i]) > max_second_delta) {
                    max_second_delta = abs(delta_second_layer[i]);
                }
                if (abs(delta_second_layer[i]) < min_second_delta) {
                    min_second_delta = abs(delta_second_layer[i]);
                }
            }

            std::vector<uint32_t> segment_index_new;
            std::vector<uint32_t> key_to_seg;
            int seg = 0;
            for (int j = 0;j < length - 2;j++) {

                segment_index_new.push_back(j);
                if (delta_second_layer[j] == 0) {
                    key_to_seg.push_back(seg);
                    while (delta_second_layer[j] == 0 && j < length - 2) {
                        key_to_seg.push_back(seg);
                        j++;

                    }

                }
                else {
                    key_to_seg.push_back(seg);
                }
                seg++;


            }
            segment_index_new.push_back(length - 2);
            key_to_seg.push_back(seg);


            for (int i = 0;i <= max_second_delta * 2;i++) {
                int aim_delta = i / 2 + 1;
                if (i % 2 == 0) { aim_delta = -aim_delta; }
                // std::cout<<"aim_delta: "<<aim_delta<<std::endl;
                for (int j = 0;j < length - 2;j++) {
                    if (delta_second_layer[j] == aim_delta) {

                        // int segment_id = key_to_seg[j];
                        int segment_id = lower_bound(j, segment_index_new.size(), segment_index_new);

                        int cost_add = 0;
                        int cost_decline = sizeof(T) * 8;
                        if (aim_delta < 0) {
                            // update delta_second_layer[end_index] of the next segment
                            int start_index = segment_index_new[segment_id + 1];
                            int end_index = segment_index_new[segment_id + 2] - 1;
                            cost_add = (end_index - start_index + 1) * abs(aim_delta);
                            if (cost_add < overhead) {
                                // merge two segments
                                if (end_index <= length) {
                                    delta_second_layer[end_index] += aim_delta;
                                }

                                // for (int idx = segment_index_new[segment_id + 1];idx < key_to_seg.size();idx++) {
                                //     key_to_seg[idx]--;
                                // }
                                if (segment_id + 1 < segment_index_new.size()) {
                                    segment_index_new.erase(segment_index_new.begin() + segment_id + 1);
                                }


                            }

                        }
                        else {
                            // update delta_second_layer[start_index-1] of the next segment
                            int start_index = segment_index_new[segment_id];
                            int end_index = segment_index_new[segment_id + 1] - 1;
                            cost_add = (end_index - start_index + 1) * abs(aim_delta);
                            if (cost_add < overhead) {
                                // merge two segments
                                if (start_index - 1 >= 0) {
                                    delta_second_layer[start_index - 1] += aim_delta;
                                }

                                // for (int idx = segment_index_new[segment_id + 1];idx < key_to_seg.size();idx++) {
                                //     key_to_seg[idx]--;
                                // }
                                if (segment_id + 1 < segment_index_new.size()) {
                                    segment_index_new.erase(segment_index_new.begin() + segment_id + 1);
                                }

                            }

                        }




                    }

                }

                // std::cout << "********" << aim_delta << "********" << std::endl;
                // for (auto item : segment_index_new) {
                //     std::cout << item << " ";
                // }
                // std::cout << std::endl;


                // for (auto item : key_to_seg) {
                //     std::cout << item << " ";
                // }
                // std::cout << std::endl;


            }


            total_byte = 0;
            int segment_total = segment_index_new.size();
            segment_index_new.push_back(nvalue * block_size + length);
            for (int i = 0;i < segment_total;i++) {
                segment_index_new[i] += nvalue * block_size;
            }
            for (int i = 0;i < segment_total;i++) {
                newsegment(segment_index_new[i], segment_index_new[i + 1] - 1);
            }

            int iter = 0;
            uint64_t cost_decline = total_byte;
            while (cost_decline > 0) {

                iter++;
                cost_decline = total_byte;
                merge(nvalue);

                double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
                // std::cout << "try "<<iter<<" segment number "<<(int)block_start_vec.size()<<" resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (sizeof(T) * block_size * 1.0);
                if (cost_decline_percent < 0.01) {
                    break;
                }

            }
            double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
            // std::cout << "segment number " << (int)block_start_vec.size() << " resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;

            for (auto item : block_start_vec) {
                block_start_vec_total.push_back(item);
            }
            for (auto item : segment_index) {
                segment_index_total.push_back(item);
            }
            for (auto item : segment_length) {
                segment_length_total.push_back(item);
            }
            total_byte_total += total_byte;
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
            segment_index.push_back((nvalue + 1) * block_size);
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
                    start_index = block_size * (nvalue + 1);
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
            segment_index.push_back(block_size * (nvalue + 1));
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
                    start_index = block_size * (nvalue + 1);
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
            new_segment_index.push_back(block_size * (nvalue + 1));
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

        T* decodeArray8(const size_t length, T* out, size_t nvalue) {
            T* res = out;
            //start_index + bit + theta0 + theta1 + numbers + delta
            segment_index_total.push_back(length);
            uint8_t maxerror;
            for (int i = 0;i < block_start_vec_total.size();i++) {
                int segment_length = segment_index_total[i + 1] - segment_index_total[i];
                uint8_t* tmpin = block_start_vec_total[i];
                tmpin += sizeof(uint32_t);
                maxerror = tmpin[0];
                tmpin++;
                if (maxerror == 127) {
                    T tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }
                if (maxerror == 126) {
                    T tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    memcpy(&tmp_val, tmpin + sizeof(T), sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }

                T base;
                memcpy(&base, tmpin, sizeof(T));
                tmpin += sizeof(T);
                if(maxerror){
                    res[0] = base;
                    read_all_bit_Delta<T>(tmpin, 0, segment_length-1, maxerror, base, res+1);
                }
                else{
                    for(int j=0;j<segment_length;j++){
                        res[j] = base;
                    }
                }

                res += segment_length;
            }
            return out;
        }



        T randomdecodeArray8(uint8_t* in, int to_find, uint32_t* out, size_t nvalue) {

            uint32_t length = segment_index_total.size();
            uint8_t* this_block = block_start_vec_total[lower_bound(to_find, length, segment_index_total)];

            uint8_t* tmpin = this_block;

            uint32_t start_ind;
            memcpy(&start_ind, tmpin, 4);
            tmpin += 4;

            T tmp_val = 0;
            uint8_t maxbits;
            memcpy(&maxbits, tmpin, sizeof(uint8_t));
            tmpin += sizeof(uint8_t);
            if (maxbits == 127) {
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                return tmp_val;
            }
            if (maxbits == 126) {
                if (to_find - start_ind == 0) {
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                }
                else {
                    tmpin += sizeof(tmp_val);
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                }
                return tmp_val;
            }

            if (maxbits == sizeof(T) * 8) {
                T tmp_val = reinterpret_cast<T*>(tmpin)[to_find];
                return tmp_val;
            }

            T base;
            memcpy(&base, tmpin, sizeof(T));
            tmpin += sizeof(T);

            tmp_val = base;
            if (maxbits) {
                tmp_val = read_Delta_int(tmpin, maxbits, to_find - start_ind, base);
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
        uint32_t get_total_byte() {
            return total_byte_total;
        }
        uint32_t get_total_blocks() {
            return block_start_vec_total.size();
        }


    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
