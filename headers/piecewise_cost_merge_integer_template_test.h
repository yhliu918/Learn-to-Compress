
#ifndef PIECEWISE_COST_INTEGER_MERGE_TEMPLATE_TEST_H_
#define PIECEWISE_COST_INTEGER_MERGE_TEMPLATE_TEST_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff
#include "stx-btree/btree.h"
#include "stx-btree/btree_map.h"
#include "ALEX/alex.h"
#include "art/art32.h"

namespace Codecset {
    template <typename T>
    class Leco_cost_merge_test
    {
    public:

        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        std::vector<uint8_t*> block_start_vec_total;
        std::vector<uint32_t> segment_index_total;
        std::vector<uint32_t> segment_length_total;

        std::vector<KeyValue<uint32_t>> art_build_vec;
        std::vector<std::pair<int, int>> alex_build_vec;
        std::vector<ART32::Node *> search_node;

        double split_time = 0;
        double merge_time = 0;
        


        uint64_t total_byte_total = 0;
        uint64_t total_byte = 0;
        int overhead = 0;
        T* array;
        int block_num;
        int block_size;
        int segment_index_total_idx = 0;

        //start_index + bit + theta0 + theta1 + numbers + delta
        void init(int blocks, int blocksize, uint64_t delta) {
            block_num = blocks;
            block_size = blocksize;
            overhead = delta; // add some punishing item

        }

        // stx::btree_map<int, int> btree_total;
        alex::Alex<int, int> alex_tree;
        ART32 art;

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



        int binarySearch2(int l, int r, int32_t x, std::vector<uint32_t>& a) {
            int d = r - l;
            int half = d >> 1;
            while (half > 0) {
                int m = l + half;
                l = (a[m] <= x) ? m : l;
                d -= half;
                half = d >> 1;
            }
            return l;
        }

        int simdSearch(int l, int r, int32_t x,  std::vector<uint32_t>& a) {
            __m256i reg2 = _mm256_set1_epi32(x);
            for (int i = l; i < r; i += 8) {
                __m256i reg1 = _mm256_load_si256((__m256i*) & a[i]);
                __m256i cmp = _mm256_cmpeq_epi32(reg1, reg2);
                unsigned bitmask = _mm256_movemask_epi8(cmp);
                if (bitmask > 0)
                    return i + (__builtin_ctz(bitmask) >> 2);
            }
            return -1;
        }

        uint64_t newsegment_size(uint32_t origin_index, uint32_t end_index) {
            
            if (origin_index == end_index) {
                return 9;
            }
            if (end_index == origin_index + 1) {
                return 13;
            }
            uint64_t overhead = sizeof(float)*2 + 5;
            int length = end_index - origin_index + 1;

            lr_int_T<T> mylr;
            mylr.caltheta(array+origin_index, length);
            float final_slope = mylr.theta1;
            float theta0 = mylr.theta0;

            int64_t max_error_delta = INT64_MIN;
            int64_t min_error_delta = INT64_MAX;
            for (int j = origin_index;j <= end_index;j++) {
                int64_t tmp = array[j] - (long long)(theta0 + final_slope * (double)(j - origin_index));
                if (tmp > max_error_delta) {
                    max_error_delta = tmp;
                }
                if (tmp < min_error_delta) {
                    min_error_delta = tmp;
                }
            }
            theta0 += (max_error_delta + min_error_delta) / 2.0;

            T final_max_error = 0;
            std::vector<bool> signvec;
            std::vector<T> delta_final;

            for (int j = origin_index;j <= end_index;j++) {
                T tmp_val;
                int128_t pred = theta0 + final_slope * (double)(j - origin_index);
                if (array[j] > pred)
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
            if (final_max_error) {
                delta_final_max_bit = bits_int_T<T>(final_max_error) + 1;
            }


            if (delta_final_max_bit >= sizeof(T) * 8) {
                delta_final_max_bit = sizeof(T) * 8;
                overhead = 5 + sizeof(T) * length;
                return overhead;
            }

            overhead += ceil((delta_final_max_bit * length)/8.0);

            return overhead;

        }

        
        void newsegment(uint32_t origin_index, uint32_t end_index) {
            
            if (origin_index == end_index) {
                return newsegment_1(origin_index, origin_index);
            }
            if (end_index == origin_index + 1) {
                return newsegment_2(origin_index, end_index);
            }
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(T) * 4);
            uint8_t* out = descriptor;
            int length = end_index - origin_index + 1;

            lr_int_T<T> mylr;
            mylr.caltheta(array+origin_index, length);
            float final_slope = mylr.theta1;
            float theta0 = mylr.theta0;

            int64_t max_error_delta = INT64_MIN;
            int64_t min_error_delta = INT64_MAX;
            for (int j = origin_index;j <= end_index;j++) {
                int64_t tmp = array[j] - (long long)(theta0 + final_slope * (double)(j - origin_index));
                if (tmp > max_error_delta) {
                    max_error_delta = tmp;
                }
                if (tmp < min_error_delta) {
                    min_error_delta = tmp;
                }
            }
            theta0 += (max_error_delta + min_error_delta) / 2.0;

            T final_max_error = 0;
            std::vector<bool> signvec;
            std::vector<T> delta_final;

            for (int j = origin_index;j <= end_index;j++) {
                // long long pred = theta0 + (float)(j - origin_index) * final_slope;
                T tmp_val;
                int128_t pred = theta0 + final_slope * (double)(j - origin_index);
                if (array[j] > pred)
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
            if (final_max_error) {
                delta_final_max_bit = bits_int_T<T>(final_max_error) + 1;
            }


            if (delta_final_max_bit >= sizeof(T) * 8) {
                delta_final_max_bit = sizeof(T) * 8;
            }

            memcpy(out, &origin_index, sizeof(origin_index));
            out += sizeof(origin_index);
            out[0] = (uint8_t)delta_final_max_bit;
            out++;

            if (delta_final_max_bit == sizeof(T) * 8) {
                for (auto i = origin_index; i <= end_index; i++)
                {
                    memcpy(out, &array[i], sizeof(T));
                    out += sizeof(T);
                }
                uint64_t segment_size = out - descriptor;
                descriptor = (uint8_t*)realloc(descriptor, segment_size);
                block_start_vec_total.push_back(descriptor);
                segment_index_total.push_back(origin_index);
                segment_length_total.push_back(segment_size);
                total_byte_total += segment_size;

                return;
            }

            memcpy(out, &theta0, sizeof(theta0));
            out += sizeof(theta0);

            memcpy(out, &final_slope, sizeof(final_slope));
            out += sizeof(final_slope);


            if (delta_final_max_bit) {
                out = write_delta_int_T(delta_final, signvec, out, delta_final_max_bit, (end_index - origin_index + 1));
            }


            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t*)realloc(descriptor, segment_size);
            block_start_vec_total.push_back(descriptor);
            segment_index_total.push_back(origin_index);
            segment_length_total.push_back(segment_size);
            total_byte_total += segment_size;

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
            block_start_vec_total.push_back(descriptor);
            segment_index_total.push_back(origin_index);
            segment_length_total.push_back(segment_size);

            total_byte_total += segment_size;
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
            block_start_vec_total.push_back(descriptor);
            segment_length_total.push_back(segment_size);
            segment_index_total.push_back(origin_index);

            total_byte_total += segment_size;
        }
        
        uint32_t cal_bits(int64_t min, int64_t max) {
            int64_t range = ceil(abs(max - min) / 2.);
            uint32_t bits = 0;
            if (range) {
                bits = bits_int_T<T>(range) + 1;
            }
            return bits;
        }

        uint8_t* encodeArray8_int(T* in, const size_t length, uint8_t* res, size_t nvalue) {
            if(nvalue == 0){
                array = in;
            }
            // array = in;
            int64_t* delta_first_layer = new int64_t[length];
            int64_t* delta_second_layer = new int64_t[length];
            std::vector<uint32_t> delta_second_bits;
            uint32_t max_second_bit = 0;
            uint32_t min_second_bit = UINT32_MAX;

            for (int i = nvalue * block_size;i < nvalue * block_size + length - 1;i++) {
                int64_t delta_first = (int64_t)in[i + 1] - (int64_t)in[i];
                delta_first_layer[i - nvalue * block_size] = delta_first;
            }

            for (int i = 0;i < length - 2;i++) {
                uint32_t delta_tmp_bit = cal_bits(delta_first_layer[i], delta_first_layer[i + 1]);
                if (delta_tmp_bit > max_second_bit) {
                    max_second_bit = delta_tmp_bit;
                }
                if (delta_tmp_bit < min_second_bit) {
                    min_second_bit = delta_tmp_bit;
                }
                delta_second_bits.push_back(delta_tmp_bit);

            }

            std::vector<uint32_t> key_to_seg;
            std::vector<int64_t> segment_max_delta;
            std::vector<int64_t> segment_min_delta;


            for (int j = 0;j <= length - 2;j++) {
                segment_index.push_back(j);
                segment_max_delta.push_back(0);
                segment_min_delta.push_back(0);
            }
            segment_index.push_back(length - 1);
            segment_max_delta.push_back(0);
            segment_min_delta.push_back(0);


            double start_timer = getNow();

            for (int aim_bit = min_second_bit; aim_bit <= max_second_bit; aim_bit++) {
                // start with smaller delta, and merge around it

                for (int j = 0;j < length - 2;j++) {
                    if (delta_second_bits[j] == aim_bit) {

                        // aim_bit is the initiate bit of the two segment, first check whether merging these two segments can reduce the cost
                        // if so, scan left & scan right to merge more segments into it
                        int segment_id = lower_bound(j, segment_index.size(), segment_index);
                        int former_index = segment_index[segment_id]; // former_index ~ start_index - 1
                        int start_index = segment_index[segment_id + 1];
                        int now_index = segment_index[segment_id + 2] - 1; // start_index ~ now_index
                        int left_bit_origin = cal_bits(segment_min_delta[segment_id], segment_max_delta[segment_id]);
                        int right_bit_origin = cal_bits(segment_min_delta[segment_id + 1], segment_max_delta[segment_id + 1]);
                        // need to compare to [seg+1]-1

                        int64_t new_max_delta = std::max(segment_max_delta[segment_id], segment_max_delta[segment_id + 1]);
                        int64_t new_min_delta = std::min(segment_min_delta[segment_id], segment_min_delta[segment_id + 1]);
                        new_max_delta = std::max(new_max_delta, delta_first_layer[segment_index[segment_id + 1] - 1]);
                        new_min_delta = std::min(new_min_delta, delta_first_layer[segment_index[segment_id + 1] - 1]);
                        int new_bit = cal_bits(new_min_delta, new_max_delta);

                        int origin_cost = (start_index - former_index) * left_bit_origin + (now_index - start_index + 1) * right_bit_origin;
                        int merged_cost = new_bit * (now_index - former_index + 1);
                        if (merged_cost - origin_cost < overhead && segment_id + 1 < segment_index.size()) {
                            // merge
                            segment_index.erase(segment_index.begin() + segment_id + 1);
                            segment_max_delta.erase(segment_max_delta.begin() + segment_id + 1);
                            segment_min_delta.erase(segment_min_delta.begin() + segment_id + 1);
                            segment_max_delta[segment_id] = new_max_delta;
                            segment_min_delta[segment_id] = new_min_delta;

                        }
                        else {
                            continue;
                        }
                        // look left & look right

                        int segment_id_search_left = segment_id - 1;
                        while (segment_id_search_left >= 0) {

                            // std::cout<<"left"<<segment_id_search_left<<std::endl;
                            int left_index = segment_index[segment_id_search_left];
                            int64_t left_max_delta = std::max(segment_max_delta[segment_id_search_left], segment_max_delta[segment_id]);
                            int64_t left_min_delta = std::min(segment_min_delta[segment_id_search_left], segment_min_delta[segment_id]);
                            // need to compare to the delta between left segment and the current segment
                            left_max_delta = std::max(left_max_delta, delta_first_layer[segment_index[segment_id] - 1]);
                            left_min_delta = std::min(left_min_delta, delta_first_layer[segment_index[segment_id] - 1]);

                            int delta_new_bit = cal_bits(left_min_delta, left_max_delta);
                            int origin_left_delta_bit = cal_bits(segment_min_delta[segment_id_search_left], segment_max_delta[segment_id_search_left]);
                            int origin_right_delta_bit = cal_bits(segment_min_delta[segment_id], segment_max_delta[segment_id]);

                            int origin_cost = (segment_index[segment_id] - left_index) * origin_left_delta_bit + (segment_index[segment_id + 1] - segment_index[segment_id]) * origin_right_delta_bit;
                            int merged_cost = delta_new_bit * (segment_index[segment_id + 1] - left_index);

                            if (merged_cost - origin_cost < overhead && segment_id < segment_index.size()) {
                                // merge

                                segment_index.erase(segment_index.begin() + segment_id);
                                segment_max_delta.erase(segment_max_delta.begin() + segment_id);
                                segment_min_delta.erase(segment_min_delta.begin() + segment_id);
                                segment_max_delta[segment_id - 1] = left_max_delta;
                                segment_min_delta[segment_id - 1] = left_min_delta;
                                segment_id_search_left--;


                            }
                            else {

                                break;
                            }

                        }


                        segment_id = lower_bound(j, segment_index.size(), segment_index);
                        int segment_id_search_right = segment_id + 1;
                        now_index = segment_index[segment_id];

                        while (segment_id_search_right + 1 < segment_index.size()) {
                            // std::cout<<"right"<<segment_id_search_right<<" / "<<segment_index_new.size()<<std::endl;
                            int right_index = segment_index[segment_id_search_right + 1];
                            int64_t right_max_delta = std::max(segment_max_delta[segment_id_search_right], segment_max_delta[segment_id]);
                            int64_t right_min_delta = std::min(segment_min_delta[segment_id_search_right], segment_min_delta[segment_id]);
                            // need to compare to the delta between seg & seg+1
                            right_max_delta = std::max(right_max_delta, delta_first_layer[segment_index[segment_id_search_right] - 1]);
                            right_min_delta = std::min(right_min_delta, delta_first_layer[segment_index[segment_id_search_right] - 1]);

                            int delta_new_bit = cal_bits(right_min_delta, right_max_delta);
                            int origin_left_delta_bit = cal_bits(segment_min_delta[segment_id], segment_max_delta[segment_id]);
                            int origin_right_delta_bit = cal_bits(segment_min_delta[segment_id_search_right], segment_max_delta[segment_id_search_right]);

                            int origin_cost = (right_index - segment_index[segment_id_search_right]) * origin_right_delta_bit + (segment_index[segment_id_search_right] - now_index) * origin_left_delta_bit;
                            int merged_cost = delta_new_bit * (right_index - now_index);

                            if (merged_cost - origin_cost < overhead) {
                                // merge
                                if (segment_id + 1 < segment_index.size()) {
                                    segment_index.erase(segment_index.begin() + segment_id + 1);
                                    segment_max_delta.erase(segment_max_delta.begin() + segment_id + 1);
                                    segment_min_delta.erase(segment_min_delta.begin() + segment_id + 1);
                                    segment_max_delta[segment_id] = right_max_delta;
                                    segment_min_delta[segment_id] = right_min_delta;
                                }
                            }
                            else {

                                break;
                            }

                        }
                        segment_id = lower_bound(j, segment_index.size(), segment_index);
                        j = segment_index[segment_id + 1] - 1;




                    }

                }

            }

            double end_timer = getNow();
            split_time +=(end_timer - start_timer);

            total_byte = 0;
            int segment_total = segment_index.size();
            segment_index.push_back(nvalue * block_size + length);
            for (int i = 0;i < segment_total;i++) {
                segment_index[i] += nvalue * block_size;
            }
            for (int i = 0;i < segment_total;i++) {
                uint64_t tmp_size = newsegment_size(segment_index[i], segment_index[i + 1] - 1);
                total_byte += tmp_size;
                segment_length.emplace_back(tmp_size);
            }
            // std::cout<< total_byte <<std::endl;

            start_timer = getNow();
            int iter = 0;
            uint64_t cost_decline = total_byte;
            while (cost_decline > 0) {

                iter++;
                cost_decline = total_byte;
                // merge(nvalue);
                merge_both_direction(nvalue, length);
                double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
                // std::cout << "try " << iter << " segment number " << (int)segment_index.size() << " resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (sizeof(T) * block_size * 1.0);
                if (cost_decline_percent < 0.01) {
                    break;
                }

            }

            int segment_number = (int)segment_index.size();
            segment_index.push_back(block_size * nvalue + length);
            for (int i = 0;i < segment_number;i++) {
                // std::cout<<segment_index[i]<<std::endl;
                newsegment(segment_index[i], segment_index[i + 1] - 1);
            }
            segment_index.pop_back();

            end_timer = getNow();
            merge_time +=(end_timer - start_timer);

            for (auto item : segment_index) {
                art_build_vec.push_back((KeyValue<uint32_t>){item, segment_index_total_idx});
                // std::cout<<item<<std::endl;
                // auto tmp = btree_total.insert(std::make_pair(item, segment_index_total_idx));
                // auto tmp = alex_tree.insert(item, segment_index_total_idx);
                alex_build_vec.push_back(std::make_pair(item, segment_index_total_idx));
                segment_index_total_idx++;
            }

            if(nvalue == block_num - 1){
                // auto tmp = btree_total.insert(std::make_pair(block_num * block_size, segment_index_total_idx));
                // std::cout<<block_num * block_size<<" "<<segment_index_total_idx<<std::endl;
                // auto tmp = alex_tree.insert(block_num * block_size, segment_index_total_idx);
                art_build_vec.push_back((KeyValue<uint32_t>){block_num * block_size, segment_index_total_idx});
                alex_build_vec.push_back(std::make_pair( block_num * block_size, segment_index_total_idx));

            }
            segment_index.clear();
            segment_length.clear();
            total_byte = 0;
            delete [] delta_first_layer;
            delete [] delta_second_layer;

            return res;

        }

        void merge_both_direction(int nvalue, int length) {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = segment_index[0];  // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int total_segments = segment_index.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;

            segment_index.push_back(block_size * nvalue + length);
            std::vector<uint32_t> new_segment_index;
            std::vector<uint32_t> new_segment_length;
            new_segment_index.emplace_back(segment_index[segment_num]);
            new_segment_length.emplace_back(segment_length[segment_num]);
            totalbyte_after_merge += segment_length[segment_num];
            segment_num++;

            while (segment_num < total_segments) {
                // std::cout<<"segment_num: "<<segment_num <<" / "<<total_segments<<std::endl;

                if (segment_num == total_segments - 1) {
                    // only can try merging with former one
                    int last_merged_segment = new_segment_length.size() - 1;
                    uint32_t init_cost_front = segment_length[segment_num] + new_segment_length[last_merged_segment];
                    uint32_t merge_cost_front = newsegment_size(new_segment_index[last_merged_segment], segment_index[segment_num + 1] - 1);
                    int saved_cost_front = init_cost_front - merge_cost_front;
                    if (saved_cost_front > 0) {
                        totalbyte_after_merge -= new_segment_length[new_segment_length.size() - 1];
                        new_segment_length[new_segment_length.size() - 1] = merge_cost_front;
                        totalbyte_after_merge += merge_cost_front;
                        segment_num++;
                    }
                    else {
                        new_segment_index.emplace_back(segment_index[segment_num]);
                        new_segment_length.emplace_back(segment_length[segment_num]);
                        totalbyte_after_merge += segment_length[segment_num];
                        segment_num++;
                    }
                    break;
                }
                int last_merged_segment = new_segment_length.size() - 1;

                uint32_t init_cost_front = segment_length[segment_num] + new_segment_length[last_merged_segment];
                uint32_t merge_cost_front = newsegment_size(new_segment_index[last_merged_segment], segment_index[segment_num + 1] - 1);
                int saved_cost_front = init_cost_front - merge_cost_front;
                // std::cout<<init_cost_front<<" "<<merge_cost_front<<std::endl;

                uint32_t init_cost_back = segment_length[segment_num] + segment_length[segment_num + 1];
                uint32_t merge_cost_back = newsegment_size(segment_index[segment_num], segment_index[segment_num + 2] - 1);
                int saved_cost_back = init_cost_back - merge_cost_back;
                // std::cout<<init_cost_back<<" "<<merge_cost_back<<std::endl;

                int saved_cost = std::max(saved_cost_front, saved_cost_back);
                if (saved_cost <= 0) {
                    // do not merge
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = segment_index[segment_num + 1];
                    segment_num++;

                    continue;
                }
                if (saved_cost_back > saved_cost_front) {
                    // merge with back
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(merge_cost_back);
                    totalbyte_after_merge += merge_cost_back;
                    start_index = segment_index[segment_num + 2];
                    segment_num += 2;
                    // std::cout<<segment_num<<std::endl;

                }
                else {
                    // merge with front
                    totalbyte_after_merge -= new_segment_length[new_segment_length.size() - 1];
                    new_segment_length[new_segment_length.size() - 1] = merge_cost_front;
                    totalbyte_after_merge += merge_cost_front;
                    start_index = segment_index[segment_num + 1];
                    segment_num += 1;

                }

            }
            total_byte = totalbyte_after_merge;
            segment_index.clear();
            segment_length.clear();

            segment_index.swap(new_segment_index);
            segment_length.swap(new_segment_length);
            // std::cout<<totalbyte_after_merge<<std::endl;

        }

        T* decodeArray8(const size_t length, T* out, size_t nvalue) {
            T* res = out;
            //start_index + bit + theta0 + theta1 + numbers + delta
            segment_index_total.push_back(length);
            float theta0;
            float theta1;
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


                memcpy(&theta0, tmpin, sizeof(theta0));
                tmpin += sizeof(theta0);
                memcpy(&theta1, tmpin, sizeof(theta1));
                tmpin += sizeof(theta1);
                if (maxerror) {
                    if (maxerror >= sizeof(T) * 8 - 1) {
                        // read_all_default(tmpin, 0, 0, segment_length, maxerror, theta1, theta0, res);
                    }
                    else {
                        read_all_bit_fix<T>(tmpin, 0, 0, segment_length, maxerror, theta1, theta0, res);
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

        int get_segment_id(int to_find) {
            // int segment_id = art.upper_bound_new(to_find, search_node) - 1;
            // int segment_id = lower_bound(to_find, segment_index_total.size(), segment_index_total);
            int segment_id = alex_tree.upper_bound(to_find).payload() -1;
            __builtin_prefetch(block_start_vec_total.data()+segment_id, 0, 3);
            return segment_id;
        }
        
        T randomdecodeArray8(int segment_id, uint8_t* in, int to_find, uint32_t* out, size_t nvalue) {


            // uint32_t length = segment_index_total.size();

            // use btree to find the segment
            // auto it = btree_total.upper_bound(to_find);
            // int segment_num = it.data();
            // uint8_t* this_block = block_start_vec_total[segment_num-1];

            // use ALEX
            // auto it = alex_tree.upper_bound(to_find);
            //  segment_id = it.payload() - 1;

            // use ART
            // s td::cout<<to_find<<std::endl;

            // int segment_id = art.upper_bound_new(to_find, search_node) - 1;

            // std::cout<<to_find<<" "<<segment_id<<std::endl;

            // normal binary search
            // segment_id = lower_bound(to_find, length, segment_index_total);
            // int segment_id = binarySearch2(0, length-1, to_find, segment_index_total);


            uint8_t* this_block = block_start_vec_total[segment_id];

            uint8_t* tmpin = this_block;

            uint32_t start_ind;
            memcpy(&start_ind, tmpin, 4);
            tmpin += 4;

            uint8_t maxerror;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            if (maxerror == sizeof(T) * 8) {
                T tmp_val = reinterpret_cast<T*>(tmpin)[to_find];
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

            float theta0;
            memcpy(&theta0, tmpin, sizeof(theta0));
            tmpin += sizeof(theta0);

            float theta1;
            memcpy(&theta1, tmpin, sizeof(theta1));
            tmpin += sizeof(theta1);


            if (maxerror) {
                // tmp_val = read_bit_fix_int_float<T>(tmpin, maxerror, to_find-start_ind, theta1, theta0);
                // tmp_val = read_bit_fix_int_float<T>(tmpin, maxerror, to_find - start_ind, theta1, theta0);
                tmp_val = read_bit_fix_T(tmpin, maxerror, to_find - start_ind, (double)theta1, theta0, 0);
            }
            else {
                tmp_val = (T)(theta0 + theta1 * (double)(to_find - start_ind));
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
            std::cout<<"split time "<<split_time<<std::endl;
            std::cout<<"merge time "<<merge_time<<std::endl;
            return block_start_vec_total.size();
        }


    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
