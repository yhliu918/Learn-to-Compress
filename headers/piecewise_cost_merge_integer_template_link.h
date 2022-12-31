
#ifndef PIECEWISE_COST_INTEGER_MERGE_TEMPLATE_TEST_LINK_H_
#define PIECEWISE_COST_INTEGER_MERGE_TEMPLATE_TEST_LINK_H_

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
    template <typename S>
    struct Segment {
        // [start, end], this maintains delta information
        int start_index;
        int end_index;
        S max_delta;
        S min_delta;
        S next_delta; // this is not contained in the segment
        int double_delta_next;
        Segment* prev;
        Segment* next;

        Segment(int start, int end, S max, S min, S next, int bit_next) {
            start_index = start;
            end_index = end;
            max_delta = max;
            min_delta = min;
            next_delta = next;
            double_delta_next = bit_next;
        }
    };

    template <typename T>
    class Leco_cost_merge_test_link
    {
    public:

        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        std::vector<uint8_t*> block_start_vec_total;
        std::vector<uint32_t> segment_index_total;
        std::vector<uint32_t> segment_length_total;

        std::vector<KeyValue<uint32_t>> art_build_vec;
        std::vector<std::pair<int, int>> alex_build_vec;
        std::vector<ART32::Node*> search_node;

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


        uint64_t newsegment_size(uint32_t origin_index, uint32_t end_index) {

            if (origin_index == end_index) {
                return 9;
            }
            if (end_index == origin_index + 1) {
                return 13;
            }
            uint64_t overhead = sizeof(float) * 2 + 5;
            int length = end_index - origin_index + 1;

            lr_int_T<T> mylr;
            mylr.caltheta(array + origin_index, length);
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

            overhead += ceil((delta_final_max_bit * length) / 8.0);

            return overhead;

        }


        void newsegment(uint32_t origin_index, uint32_t end_index) {

            // if(origin_index<=127452 && 127452<=end_index){
            //     std::cout<<"hello"<<std::endl;
            // }
            if (origin_index == end_index) {
                return newsegment_1(origin_index, origin_index);
            }
            if (end_index == origin_index + 1) {
                return newsegment_2(origin_index, end_index);
            }
            // if(origin_index==9980){
            //     std::cout<<"heloo"<<std::endl;
            // }
            // std::cout<<origin_index<<" "<<(end_index - origin_index + 1) <<std::endl;
            uint8_t* descriptor = (uint8_t*)malloc((end_index - origin_index + 1) * sizeof(T) * 4+200);
            uint8_t* out = descriptor;
            int length = end_index - origin_index + 1;

            lr_int_T<T> mylr;
            mylr.caltheta(array + origin_index, length);
            float final_slope = mylr.theta1;
            double theta0 = mylr.theta0;

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
            if(abs(final_slope) >= 0.00000001 && delta_final_max_bit != sizeof(T) * 8){ 
                out[0] += (1<<7); // if first bit of out[0] is 1 means have slope, else slope = 0
            }
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
            if(abs(final_slope) >= 0.00000001){ 
                memcpy(out, &final_slope, sizeof(final_slope));
                out += sizeof(final_slope);
            }

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
            out[0] = (uint8_t)254; // this means that this segment only has two points
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
            out[0] = (uint8_t)255; // this means that this segment only has one point
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
            array = in;

            Segment<int64_t> head(0, 0, 0, 0, 0, 10000);
            Segment<int64_t> tail(0, 0, 0, 0, 0, 0);
            Segment<int64_t>* current = &head;
            int min_second_bit = 10000;
            int max_second_bit = -1;
            int64_t delta_prev = int64_t(array[nvalue * block_size + 1]) - int64_t(array[nvalue * block_size]);
            for (int i = nvalue * block_size + 1;i < nvalue * block_size + length - 1;i++) {
                int64_t delta = int64_t(array[i + 1]) - int64_t(array[i]);
                int second_delta_bit = cal_bits(delta_prev, delta);
                if (second_delta_bit < min_second_bit) {
                    min_second_bit = second_delta_bit;
                }
                if (second_delta_bit > max_second_bit) {
                    max_second_bit = second_delta_bit;
                }
                Segment <int64_t>* newseg = new Segment<int64_t>(i - 1, i - 1, 0, 0, delta_prev, second_delta_bit);
                current->next = newseg;
                newseg->prev = current;
                current = newseg;
                delta_prev = delta;
            }
            Segment<int64_t>* newseg = new Segment<int64_t>(nvalue * block_size + length - 2, nvalue * block_size + length - 2, 0, 0, delta_prev, 10000);
            current->next = newseg;
            newseg->prev = current;
            current = newseg;
            current->next = &tail;
            tail.prev = current;

            // current = (&head)->next;
            // while (current->next != &tail) {
            //     if(current->delta_bit_next!=0){
            //             std::cout<<current->start_index<<" "<<current->end_index<<" "<<current->delta_bit_next<<std::endl;
            //     }

            //     current = current->next;
            // }
            bool flag = false;

            double start_timer = getNow();

            for (int aim_bit = min_second_bit; aim_bit <= max_second_bit;aim_bit++) {
                current = (&head)->next;
                while (current != &tail && current->next != &tail) {

                    if (current->double_delta_next == aim_bit) {
                        Segment<int64_t>* next = current->next;
                        int former_index = current->start_index; // former_index ~ start_index - 1
                        int start_index = next->start_index;
                        int now_index = next->end_index; // start_index ~ now_index
                        int left_bit_origin = cal_bits(current->min_delta, current->max_delta);
                        int right_bit_origin = cal_bits(next->min_delta, next->max_delta);

                        int64_t new_max_delta = std::max(current->max_delta, next->max_delta);
                        new_max_delta = std::max(new_max_delta, current->next_delta);
                        int64_t new_min_delta = std::min(current->min_delta, next->min_delta);
                        new_min_delta = std::min(new_min_delta, current->next_delta);
                        int new_bit = cal_bits(new_min_delta, new_max_delta);

                        int origin_cost = (start_index - former_index) * left_bit_origin + (now_index - start_index + 1) * right_bit_origin;
                        int merged_cost = new_bit * (now_index - former_index + 1);
                        if (merged_cost - origin_cost < overhead) {
                            // merge
                            current->end_index = now_index;
                            current->next = next->next;
                            next->next->prev = current;
                            current->next_delta = next->next_delta;
                            current->max_delta = new_max_delta;
                            current->min_delta = new_min_delta;
                            current->double_delta_next = next->double_delta_next;
                            delete next;

                        }
                        else {
                            current = current->next;
                            continue;
                        }
                        // look left
                        Segment<int64_t>* prev = current->prev;
                        while (prev != &head && prev->prev != &head) {
                            int left_index = prev->start_index;
                            int64_t left_max_delta = std::max(prev->max_delta, current->max_delta);
                            left_max_delta = std::max(left_max_delta, prev->next_delta);
                            int64_t left_min_delta = std::min(prev->min_delta, current->min_delta);
                            left_min_delta = std::min(left_min_delta, prev->next_delta);

                            int new_bit = cal_bits(left_min_delta, left_max_delta);
                            int origin_left_delta_bit = cal_bits(prev->min_delta, prev->max_delta);
                            int origin_right_delta_bit = cal_bits(current->min_delta, current->max_delta);
                            int origin_cost = (current->start_index - left_index) * origin_left_delta_bit + (current->end_index - current->start_index + 1) * origin_right_delta_bit;
                            int merged_cost = new_bit * (current->end_index - left_index + 1);

                            if (merged_cost - origin_cost < overhead) {
                                // merge
                                current->start_index = left_index;
                                current->prev = prev->prev;
                                prev->prev->next = current;
                                current->min_delta = left_min_delta;
                                current->max_delta = left_max_delta;
                                delete prev;
                                prev = current->prev;
                            }
                            else {

                                break;
                            }

                        }

                        next = current->next;
                        while (next != &tail && next->next != &tail) {
                            int right_index = next->end_index;
                            int64_t right_max_delta = std::max(next->max_delta, current->max_delta);
                            right_max_delta = std::max(right_max_delta, current->next_delta);
                            int64_t right_min_delta = std::min(next->min_delta, current->min_delta);
                            right_min_delta = std::min(right_min_delta, current->next_delta);

                            int new_bit = cal_bits(right_min_delta, right_max_delta);
                            int origin_left_delta_bit = cal_bits(current->min_delta, current->max_delta);
                            int origin_right_delta_bit = cal_bits(next->min_delta, next->max_delta);
                            int origin_cost = (right_index - next->start_index + 1) * origin_right_delta_bit + (next->start_index - current->start_index) * origin_left_delta_bit;
                            int merged_cost = new_bit * (right_index - current->start_index + 1);

                            if (merged_cost - origin_cost < overhead) {
                                // merge
                                current->end_index = right_index;
                                current->next = next->next;
                                next->next->prev = current;
                                current->max_delta = right_max_delta;
                                current->min_delta = right_min_delta;
                                current->double_delta_next = next->double_delta_next;
                                current->next_delta = next->next_delta;
                                delete next;
                                next = current->next;
                            }
                            else {
                                break;
                            }

                        }

                        current = current->next;
                    }
                    else{
                        current = current->next;
                    }

                }
            }


            double end_timer = getNow();
            split_time += (end_timer - start_timer);

            current = (&head)->next;
            while (current->next != &tail) {
                segment_index.push_back(current->start_index );
                current = current->next;
            }

            int segment_total = segment_index.size();
            segment_index.push_back(nvalue * block_size + length);

            total_byte = 0;
            for (int i = 0;i < segment_total;i++) {
                uint64_t tmp_size = newsegment_size(segment_index[i], segment_index[i + 1] - 1);
                total_byte += tmp_size;
                segment_length.emplace_back(tmp_size);
            }
            segment_index.pop_back();
            // std::cout << total_byte << std::endl;

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
            if(segment_number<=1){
                segment_index.clear();
                segment_index.push_back(block_size * nvalue);
                segment_index.push_back(block_size * nvalue + length);
                newsegment(block_size * nvalue, block_size * nvalue+length - 1);
            }
            else{
                segment_index.push_back(block_size * nvalue + length);
                for (int i = 0;i < segment_number;i++) {
                // std::cout<<segment_index[i]<<std::endl;
                newsegment(segment_index[i], segment_index[i + 1] - 1);
                }
            }
            
            segment_index.pop_back();

            end_timer = getNow();
            merge_time += (end_timer - start_timer);

            for (auto item : segment_index) {
                art_build_vec.push_back((KeyValue<uint32_t>) { item, segment_index_total_idx });
                // std::cout<<item<<std::endl;
                // auto tmp = btree_total.insert(std::make_pair(item, segment_index_total_idx));
                // auto tmp = alex_tree.insert(item, segment_index_total_idx);
                alex_build_vec.push_back(std::make_pair(item, segment_index_total_idx));
                segment_index_total_idx++;
            }

            if (nvalue == block_num - 1) {
                // auto tmp = btree_total.insert(std::make_pair(block_num * block_size, segment_index_total_idx));
                // std::cout<<block_num * block_size<<" "<<segment_index_total_idx<<std::endl;
                // auto tmp = alex_tree.insert(block_num * block_size, segment_index_total_idx);
                art_build_vec.push_back((KeyValue<uint32_t>) { block_num* block_size, segment_index_total_idx });
                alex_build_vec.push_back(std::make_pair(block_num * block_size, segment_index_total_idx));

            }
            segment_index.clear();
            segment_length.clear();
            total_byte = 0;
            current = (&head)->next;
            Segment<int64_t>* next = current->next;
            while (current != &tail && current->next != &tail) {
                next = current->next;
                delete current;
                current = next;
            }
            (&head)->next = &tail;
            (&tail)->prev = &head;

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
            double theta0 = 0;
            float theta1 = 0;
            uint8_t maxerror;
            for (int i = 0;i < block_start_vec_total.size();i++) {
                int segment_length = segment_index_total[i + 1] - segment_index_total[i];
                uint8_t* tmpin = block_start_vec_total[i];
                //debug
                int start_ind = 0;
                memcpy(&start_ind, tmpin, sizeof(int));

                //debug
                tmpin += sizeof(uint32_t);
                maxerror = tmpin[0];
                tmpin++;
                if (maxerror == 255) {
                    T tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }
                if (maxerror == 254) {
                    T tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    memcpy(&tmp_val, tmpin + sizeof(T), sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }
                if(maxerror==sizeof(T)*8){
                    memcpy(res, tmpin, sizeof(T)*segment_length);
                    res+=segment_length;
                    continue;
                }
                
                memcpy(&theta0, tmpin, sizeof(theta0));
                tmpin += sizeof(theta0);
                theta1 = 0;
                if((maxerror>>7) == 1){
                    memcpy(&theta1, tmpin, sizeof(theta1));
                    tmpin += sizeof(theta1);
                    maxerror -= 128;
                } 
                if (maxerror) {
                    read_all_bit_fix<T>(tmpin, 0, 0, segment_length, maxerror, theta1, theta0, res);
                    
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
            int segment_id = alex_tree.upper_bound(to_find).payload() - 1;
            __builtin_prefetch(block_start_vec_total.data() + segment_id, 0, 3);
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
                T tmp_val = reinterpret_cast<T*>(tmpin)[to_find-start_ind];
                return tmp_val;
            }

            T tmp_val = 0;
            if (maxerror == 255) {
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                return tmp_val;
            }
            if (maxerror == 254) {
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

            float theta1=0;
            if((maxerror>>7) == 1){
                memcpy(&theta1, tmpin, sizeof(theta1));
                tmpin += sizeof(theta1);
                maxerror -= 128;
            } 


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
            // std::cout << "split time " << split_time << std::endl;
            // std::cout << "merge time " << merge_time << std::endl;
            return block_start_vec_total.size();
        }


    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
