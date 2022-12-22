
#ifndef DELTA_COST_INTEGER_MERGE_TEMPLATE_LINK_H_
#define DELTA_COST_INTEGER_MERGE_TEMPLATE_LINK_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {
    template <typename S>
    struct Segment {
        // [start, end], this maintains delta information
        int start_index;
        int end_index;
        int delta_bit;
        int next_delta_bit; // this is not contained in the segment
        int double_delta_next;
        Segment* prev;
        Segment* next;

        Segment(int start, int end, int now, int next, int bit_next) {
            start_index = start;
            end_index = end;
            delta_bit = now;
            next_delta_bit = next;
            double_delta_next = bit_next;
        }
    };

    template <typename T>
    class Delta_cost_merge_link
    {
    public:

        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        std::vector<uint8_t*> block_start_vec_total;
        std::vector<uint32_t> segment_index_total;
        std::vector<uint32_t> segment_length_total;

        double split_time = 0;
        double merge_time = 0;

        uint64_t total_byte_total = 0;
        uint64_t total_byte = 0;
        int tau = 0;
        T* array;
        int block_num;
        int block_size;
        int segment_index_total_idx = 0;

        //start_index + bit + theta0 + theta1 + numbers + delta
        void init(int blocks, int blocksize, uint64_t delta) {
            block_num = blocks;
            block_size = blocksize;
            tau = delta; // add some punishing item

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


        uint64_t newsegment_size(uint32_t origin_index, uint32_t end_index) {

            if (origin_index == end_index) {
                return 5+sizeof(T);
            }
            if (end_index == origin_index + 1) {
                return 5+2*sizeof(T);
            }
            uint64_t overhead = sizeof(T) + 5;
            int length = end_index - origin_index + 1;

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
                overhead = 5 + sizeof(T) * length;
                return overhead;
            }

            overhead += ceil((max_bit * length) / 8.0);

            return overhead;

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
                // std::cout<<origin_index<<" "<<end_index<<" "<<unsigned(max_bit)<<std::endl;
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
                block_start_vec_total.push_back(descriptor);
                segment_index_total.push_back(origin_index);
                segment_length_total.push_back(segment_size);
                total_byte_total += segment_size;
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
            block_start_vec_total.push_back(descriptor);
            segment_index_total.push_back(origin_index);
            segment_length_total.push_back(segment_size);
            total_byte_total += segment_size;
            // std::cout<<origin_index<<" "<<total_byte_total<<std::endl;
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
            segment_index_total.push_back(origin_index);
            segment_length_total.push_back(segment_size);
            total_byte_total += segment_size;
        }

        uint8_t* encodeArray8_int(T* in, const size_t length, uint8_t* res, size_t nvalue) {
            array = in;

            Segment<int64_t> head(0, 0, 0, 0, 10000);
            Segment<int64_t> tail(0, 0, 0, 0, 0);
            Segment<int64_t>* current = &head;
            int min_second_bit = 10000;
            int max_second_bit = -1;
            int64_t delta_prev_bit = bits_int_T(array[nvalue * block_size + 1] - array[nvalue * block_size]);
            for (int i = nvalue * block_size + 1;i < nvalue * block_size + length - 1;i++) {
                int64_t delta_bit = 0;
                if(array[i+1]<array[i]){
                    delta_bit = bits_int_T(array[i] - array[i + 1]);
                }
                else{
                    delta_bit = bits_int_T(array[i + 1] - array[i]);
                }
                int second_delta_bit = delta_bit - delta_prev_bit;
                if (second_delta_bit < min_second_bit) {
                    min_second_bit = second_delta_bit;
                }
                if (second_delta_bit > max_second_bit) {
                    max_second_bit = second_delta_bit;
                }
                Segment <int64_t>* newseg = new Segment<int64_t>(i - 1, i - 1, 0, delta_prev_bit, second_delta_bit);
                current->next = newseg;
                newseg->prev = current;
                current = newseg;
                delta_prev_bit = delta_bit;
            }
            Segment<int64_t>* newseg = new Segment<int64_t>(nvalue * block_size + length - 2, nvalue * block_size + length - 2, 0, delta_prev_bit, 10000);
            current->next = newseg;
            newseg->prev = current;
            current = newseg;
            current->next = &tail;
            tail.prev = current;
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
                        int left_bit_origin = current->delta_bit;
                        int right_bit_origin = next->delta_bit;

                        int new_max_delta = std::max(current->delta_bit, next->delta_bit);
                        new_max_delta = std::max(new_max_delta, current->next_delta_bit);

                        int origin_cost = (start_index - former_index) * left_bit_origin + (now_index - start_index + 1) * right_bit_origin;
                        int merged_cost = new_max_delta * (now_index - former_index + 1);
                        if (merged_cost - origin_cost < tau) {
                            // merge
                            current->end_index = now_index;
                            current->next = next->next;
                            next->next->prev = current;
                            current->delta_bit = new_max_delta;
                            current->next_delta_bit = next->next_delta_bit;
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
                            int left_max_delta = std::max(prev->delta_bit, current->delta_bit);
                            left_max_delta = std::max(left_max_delta, prev->next_delta_bit);

                            int origin_left_delta_bit = prev->delta_bit;
                            int origin_right_delta_bit = current->delta_bit;
                            int origin_cost = (current->start_index - left_index) * origin_left_delta_bit + (current->end_index - current->start_index + 1) * origin_right_delta_bit;
                            int merged_cost = left_max_delta * (current->end_index - left_index + 1);

                            if (merged_cost - origin_cost < tau) {
                                // merge
                                current->start_index = left_index;
                                current->prev = prev->prev;
                                prev->prev->next = current;
                                current->delta_bit = left_max_delta;
                                delete prev;
                                prev = current->prev;
                            }
                            else {

                                break;
                            }

                        }
                        // look right
                        next = current->next;
                        while (next != &tail && next->next != &tail) {
                            int right_index = next->end_index;
                            int right_max_delta = std::max(current->delta_bit, next->delta_bit);
                            right_max_delta = std::max(right_max_delta, current->next_delta_bit);
                            int origin_left_delta_bit = current->delta_bit;
                            int origin_right_delta_bit = next->delta_bit;
                            int origin_cost = (right_index - next->start_index + 1) * origin_right_delta_bit + (next->start_index - current->start_index) * origin_left_delta_bit;
                            int merged_cost = right_max_delta * (right_index - current->start_index + 1);

                            if (merged_cost - origin_cost < tau) {
                                // merge
                                current->end_index = right_index;
                                current->next = next->next;
                                next->next->prev = current;
                                current->delta_bit = right_max_delta;
                                current->double_delta_next = next->double_delta_next;
                                current->next_delta_bit = next->next_delta_bit;
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
                // std::cout<<current->start_index<<std::endl;
                segment_index.push_back(current->start_index );
                current = current->next;
            }
            if(current->next == &tail){
                // std::cout<<current->start_index<<std::endl;
                segment_index.push_back(current->start_index);
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
                merge(nvalue, length);
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
            merge_time += (end_timer - start_timer);


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


        void merge(int nvalue, int length) {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = segment_index[0]; // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int newsegment_num = 0;
            int total_segments = segment_index.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;
            segment_index.push_back(nvalue  * block_size + length);
     
            std::vector<uint32_t> new_segment_index;
            std::vector<uint32_t> new_segment_length;
            while (segment_num < total_segments) {
                // std::cout<<"segment_num: "<<segment_num <<" / "<<total_segments<<std::endl;

                if (segment_num == total_segments - 1) {
                    // std::cout <<segment_num<<"///"<<total_segments<<" "<< block_start_vec[segment_num] << std::endl;
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = block_size * (nvalue + 1);
                    segment_num++;
                    break;
                }
                uint32_t init_cost = segment_length[segment_num] + segment_length[segment_num + 1];
                uint32_t merge_cost = newsegment_size(start_index, segment_index[segment_num + 2] - 1);
                if (init_cost > merge_cost) { // merge the two segments
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
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = segment_index[segment_num + 1];
                    segment_num++;
                    newsegment_num++;
                }

            }
            segment_index.swap(new_segment_index);
            segment_length.swap(new_segment_length);
            total_byte = totalbyte_after_merge;
            // std::cout<<total_byte<<std::endl;

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
                if(maxerror==sizeof(T)*8){
                    memcpy(res, tmpin, sizeof(T)*segment_length);
                    res+=segment_length;
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
                T tmp_val = reinterpret_cast<T*>(tmpin)[to_find-start_ind];
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
            // std::cout << "split time " << split_time << std::endl;
            // std::cout << "merge time " << merge_time << std::endl;
            return block_start_vec_total.size();
        }


    };



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
