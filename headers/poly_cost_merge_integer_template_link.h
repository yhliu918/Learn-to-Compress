
#ifndef POLY_COST_INTEGER_MERGE_TEMPLATE_TEST_LINK_H_
#define POLY_COST_INTEGER_MERGE_TEMPLATE_TEST_LINK_H_

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
#include "sdlp.hpp"
#include "polynomial_regression.hpp"
using namespace Eigen;

namespace Codecset
{
    template <typename S>
    struct Segment_pol
    {
        // [start, end], this maintains delta information
        int start_index;
        int end_index;
        int double_delta_next;
        int estimate_bit;
        Segment_pol *prev;
        Segment_pol *next;

        Segment_pol(int start, int end, int bit_next, int est_bit)
        {
            start_index = start;
            end_index = end;
            double_delta_next = bit_next;
            estimate_bit = est_bit;
        }
    };
    template <typename T, int degree>
    int64_t solve_lp(T* y_val, int number)
    {
        if(number<=degree+1){
            return 0;
        }
        int m = number * 2;
        Eigen::Matrix<double, degree + 2, 1> x;                    
        Eigen::Matrix<double, degree + 2, 1> c;                   
        Eigen::Matrix<double, -1, degree + 2> A(m + 1, degree + 2);
        Eigen::VectorXd b(m + 1);                                
        for (int i = 0; i < degree + 1; i++)
        {
            c(i) = 0;
        }
        c(degree + 1) = 1.0;
        for (int i = 0; i < number; i++)
        {
            for (int j = 0; j < degree + 1; j++)
            {
                A(2 * i, j) = pow(i, degree - j);
                A(2 * i + 1, j) = -pow(i, degree - j);
            }
            A(2 * i, degree + 1) = -1;
            A(2 * i + 1, degree + 1) = -1;
            b(2 * i) = (int64_t)y_val[i];
            b(2 * i + 1) = -(int64_t)y_val[i];
        }
        for (int i = 0; i < degree + 1; i++)
        {
            A(m, i) = 0;
        }
        A(m, degree + 1) = -1;
        b(m) = 0;
        double minobj = sdlp::linprog<degree + 2>(c, A, b, x);
        return int64_t(ceil(minobj));
    }

    template <typename T, int degree>
    class Poly_cost_merge_test_link
    {
    public:
        std::vector<uint32_t> segment_index;
        std::vector<uint32_t> segment_length;

        std::vector<uint8_t *> block_start_vec_total;
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
        T *array;
        int block_num;
        int block_size;
        int segment_index_total_idx = 0;

        // start_index + bit + theta0 + theta1 + numbers + delta
        void init(int blocks, int blocksize, uint64_t delta)
        {
            block_num = blocks;
            block_size = blocksize;
            overhead = delta; // add some punishing item
        }

        // stx::btree_map<int, int> btree_total;
        alex::Alex<int, int> alex_tree;
        ART32 art;

        uint32_t lower_bound(uint64_t v, uint32_t len, std::vector<uint32_t> &index)
        {
            uint32_t m;
            uint32_t x = 0;
            uint32_t y = len - 1;
            while (x <= y)
            {

                m = x + (y - x) / 2;
                if (v < index[m])
                    y = m - 1;
                else
                    x = m + 1;
            }
            return y;
        }

        uint64_t newsegment_size(uint32_t origin_index, uint32_t end_index)
        {
            if (origin_index == end_index)
            {
                return 9;
            }
            if (end_index == origin_index + 1)
            {
                return 13;
            }
            uint8_t *descriptor = (uint8_t *)malloc((end_index - origin_index + 1) * sizeof(T) * 4 + 200);
            uint8_t *out = descriptor;
            int length = end_index - origin_index + 1;

            std::vector<T> data_vec = std::vector<T>(array+origin_index, array + end_index + 1);
            auto simple_fixed = andviane::polynomial_regression<degree>(data_vec); 
            
            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;
            for (auto i = 0; i < length; i++)
            {
                T tmp_val;
                int64_t pred = simple_fixed(i);

                if ((int64_t)data_vec[i] > pred)
                {
                    tmp_val = (int64_t)data_vec[i] - pred;
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = pred - (int64_t)data_vec[i];
                    signvec.emplace_back(false); // means negative
                }
                delta.emplace_back(tmp_val);

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }
            uint8_t delta_final_max_bit = 0;
            if (max_error)
            {
                delta_final_max_bit = bits_int_T(max_error) + 1;
            }

            if (delta_final_max_bit > sizeof(T) * 8)
            {
                delta_final_max_bit = sizeof(T) * 8;
            }
            int start_ind =origin_index;
            memcpy(out, &start_ind, sizeof(start_ind));
            out += sizeof(start_ind);
            memcpy(out, &delta_final_max_bit, sizeof(delta_final_max_bit));
            out += sizeof(delta_final_max_bit);


            if (delta_final_max_bit == sizeof(T) * 8)
            {
                for (auto i = origin_index; i <= end_index; i++)
                {
                    memcpy(out, &array[i], sizeof(T));
                    out += sizeof(T);
                }
                uint64_t segment_size = out - descriptor;
                

                return segment_size;
            }

            for(auto a: simple_fixed){
                memcpy(out, &a, sizeof(double));
                out += sizeof(double);
            }

            if (delta_final_max_bit)
            {
                out = write_delta_int_T(delta, signvec, out, delta_final_max_bit, (end_index - origin_index + 1));
            }

            uint64_t segment_size = out - descriptor;
            free(descriptor);
            return segment_size;
        }

        // uint64_t newsegment_size(uint32_t origin_index, uint32_t end_index)
        // {
        //     if (origin_index == end_index)
        //     {
        //         return 9;
        //     }
        //     if (end_index == origin_index + 1)
        //     {
        //         return 13;
        //     }
        //     uint64_t overhead = sizeof(double) * (degree+1) + 5;
        //     int length = end_index - origin_index + 1;
        //     int delta_final_max_bit = cal_bits_range(solve_lp<T, degree>(array + origin_index, length));
            
        //     if (delta_final_max_bit >= sizeof(T) * 8)
        //     {
        //         delta_final_max_bit = sizeof(T) * 8;
        //         overhead = 5 + sizeof(T) * length;
        //         return overhead;
        //     }

        //     overhead += ceil((delta_final_max_bit * length) / 8.0);

        //     return overhead;
        // }

        void newsegment(uint32_t origin_index, uint32_t end_index, int nvalue)
        {
            if (origin_index == end_index)
            {
                return newsegment_1(origin_index, origin_index,nvalue);
            }
            if (end_index == origin_index + 1)
            {
                return newsegment_2(origin_index, end_index,nvalue);
            }
            uint8_t *descriptor = (uint8_t *)malloc((end_index - origin_index + 1) * sizeof(T) * 4 + 200);
            uint8_t *out = descriptor;
            int length = end_index - origin_index + 1;

            std::vector<T> data_vec = std::vector<T>(array+origin_index, array + end_index + 1);
            auto simple_fixed = andviane::polynomial_regression<degree>(data_vec); 
            
            std::vector<T> delta;
            std::vector<bool> signvec;
            T max_error = 0;
            for (auto i = 0; i < length; i++)
            {
                T tmp_val;
                int64_t pred = simple_fixed(i);

                if ((int64_t)data_vec[i] > pred)
                {
                    tmp_val = (int64_t)data_vec[i] - pred;
                    signvec.emplace_back(true); // means positive
                }
                else
                {
                    tmp_val = pred - (int64_t)data_vec[i];
                    signvec.emplace_back(false); // means negative
                }
                delta.emplace_back(tmp_val);

                if (tmp_val > max_error)
                {
                    max_error = tmp_val;
                }
            }
            uint8_t delta_final_max_bit = 0;
            if (max_error)
            {
                delta_final_max_bit = bits_int_T(max_error) + 1;
            }

            if (delta_final_max_bit > sizeof(T) * 8)
            {
                delta_final_max_bit = sizeof(T) * 8;
            }
            int start_ind = nvalue * block_size+origin_index;
            memcpy(out, &start_ind, sizeof(start_ind));
            out += sizeof(start_ind);
            memcpy(out, &delta_final_max_bit, sizeof(delta_final_max_bit));
            out += sizeof(delta_final_max_bit);


            if (delta_final_max_bit == sizeof(T) * 8)
            {
                for (auto i = origin_index; i <= end_index; i++)
                {
                    memcpy(out, &array[i], sizeof(T));
                    out += sizeof(T);
                }
                uint64_t segment_size = out - descriptor;
                descriptor = (uint8_t *)realloc(descriptor, segment_size);
                block_start_vec_total.push_back(descriptor);
                segment_index_total.push_back(origin_index+nvalue*block_size);
                segment_length_total.push_back(segment_size);
                total_byte_total += segment_size;

                return;
            }

            for(auto a: simple_fixed){
                memcpy(out, &a, sizeof(double));
                out += sizeof(double);
            }

            if (delta_final_max_bit)
            {
                out = write_delta_int_T(delta, signvec, out, delta_final_max_bit, (end_index - origin_index + 1));
            }

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t *)realloc(descriptor, segment_size);
            block_start_vec_total.push_back(descriptor);
            segment_index_total.push_back(origin_index+nvalue*block_size);
            segment_length_total.push_back(segment_size);
            total_byte_total += segment_size;
        }

        void newsegment_2(uint32_t origin_index, uint32_t end_index,int nvalue)
        {

            uint8_t *descriptor = (uint8_t *)malloc((end_index - origin_index + 1) * sizeof(T) * 2);
            uint8_t *out = descriptor;
            int start_ind = nvalue * block_size+origin_index;
            memcpy(out, &start_ind, sizeof(start_ind));
            out += sizeof(start_ind);
            out[0] = (uint8_t)254; // this means that this segment only has two points
            out++;
            memcpy(out, &array[origin_index], sizeof(T));
            out += sizeof(T);
            memcpy(out, &(array[origin_index + 1]), sizeof(T));
            out += sizeof(T);

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t *)realloc(descriptor, segment_size);
            block_start_vec_total.push_back(descriptor);
            segment_index_total.push_back(origin_index+nvalue*block_size);
            segment_length_total.push_back(segment_size);

            total_byte_total += segment_size;
        }

        void newsegment_1(uint32_t origin_index, uint32_t end_index,int nvalue)
        {

            uint8_t *descriptor = (uint8_t *)malloc(10 * sizeof(T));
            uint8_t *out = descriptor;
            int start_ind = nvalue * block_size+origin_index;
            memcpy(out, &start_ind, sizeof(start_ind));
            out += sizeof(start_ind);
            out[0] = (uint8_t)255; // this means that this segment only has one point
            out++;
            memcpy(out, &array[origin_index], sizeof(T));
            out += sizeof(T);

            uint64_t segment_size = out - descriptor;
            descriptor = (uint8_t *)realloc(descriptor, segment_size);
            block_start_vec_total.push_back(descriptor);
            segment_length_total.push_back(segment_size);
            segment_index_total.push_back(origin_index+nvalue*block_size);

            total_byte_total += segment_size;
        }


        uint32_t cal_bits_range(int64_t data_range)
        {
            uint32_t bits = 0;
            if (data_range)
            {
                bits = bits_int_T<T>(data_range) + 1;
            }
            return bits;
        }

        uint8_t *encodeArray8_int(T *in, const size_t length, uint8_t *res, size_t nvalue)
        {
            array = in + nvalue * block_size;

            Segment_pol<int64_t> head(0, 0, 10000, 0);
            Segment_pol<int64_t> tail(0, 0, 0, 0);
            Segment_pol<int64_t> *current = &head;
            int min_second_bit = 10000;
            int max_second_bit = -1;
            std::vector<int64_t> first_order_delta;
            for(int i=0;i<length-1;i++){
                first_order_delta.push_back(int64_t(array[i+1])-int64_t(array[i]));
            }
            std::vector<int64_t> second_order_delta;
            for(int i=0;i<first_order_delta.size()-2;i++){
                second_order_delta.push_back(first_order_delta[i+1]-first_order_delta[i]);
            }
            std::vector<int64_t> third_order_delta;
            for(int i=0;i<first_order_delta.size()-3;i++){
                third_order_delta.push_back(second_order_delta[i+1]-second_order_delta[i]);
            }
            for (int i = 0; i <third_order_delta.size() - 1; i++)
            {
                int64_t delta = abs(third_order_delta[i + 1] - third_order_delta[i]);
                int second_delta_bit = cal_bits_range(delta);
                if (second_delta_bit < min_second_bit)
                {
                    min_second_bit = second_delta_bit;
                }
                if (second_delta_bit > max_second_bit)
                {
                    max_second_bit = second_delta_bit;
                }
                Segment_pol<int64_t> *newseg = new Segment_pol<int64_t>(i, i, second_delta_bit, 0);
                current->next = newseg;
                newseg->prev = current;
                current = newseg;
            }
            Segment_pol<int64_t> *newseg = new Segment_pol<int64_t>(length - 3,  length - 3, 10000, 0);
            current->next = newseg;
            newseg->prev = current;
            current = newseg;
            current->next = &tail;
            tail.prev = current;

            bool flag = false;

            double start_timer = getNow();

            for (int aim_bit = min_second_bit; aim_bit <= min_second_bit+1; aim_bit++)
            {
                current = (&head)->next;
                while (current != &tail && current->next != &tail)
                {

                    if (current->double_delta_next == aim_bit)
                    {
                        Segment_pol<int64_t> *next = current->next;
                        int former_index = current->start_index; // former_index ~ start_index - 1
                        int start_index = next->start_index;
                        int now_index = next->end_index; // start_index ~ now_index
                        int left_bit_origin = current->estimate_bit;
                        int right_bit_origin = next->estimate_bit;

                        int new_bit = cal_bits_range(solve_lp<T, degree>(array+former_index, now_index - former_index+1));

                        int origin_cost = (start_index - former_index) * left_bit_origin + (now_index - start_index + 1) * right_bit_origin;
                        int merged_cost = new_bit * (now_index - former_index + 1);
                        if (merged_cost - origin_cost < overhead)
                        {
                            // merge
                            current->end_index = now_index;
                            current->next = next->next;
                            next->next->prev = current;
                            current->double_delta_next = next->double_delta_next;
                            current->estimate_bit = new_bit;
                            delete next;
                        }
                        else
                        {
                            current = current->next;
                            continue;
                        }
                        // look left
                        Segment_pol<int64_t> *prev = current->prev;
                        while (prev != &head )
                        {
                            int left_index = prev->start_index;
                            int new_bit = cal_bits_range(solve_lp<T, degree>(array+left_index, current->end_index - left_index+1));
                            int origin_left_delta_bit = prev->estimate_bit;
                            int origin_right_delta_bit = current->estimate_bit;
                            int origin_cost = (current->start_index - left_index) * origin_left_delta_bit + (current->end_index - current->start_index + 1) * origin_right_delta_bit;
                            int merged_cost = new_bit * (current->end_index - left_index + 1);

                            if (merged_cost - origin_cost < overhead)
                            {
                                // merge
                                current->start_index = left_index;
                                current->prev = prev->prev;
                                prev->prev->next = current;
                                current->estimate_bit = new_bit;
                                delete prev;
                                prev = current->prev;
                            }
                            else
                            {

                                break;
                            }
                        }

                        next = current->next;
                        while (next != &tail && next->next != &tail)
                        {
                            int right_index = next->end_index;
                            int new_bit = cal_bits_range(solve_lp<T, degree>(array+current->start_index, right_index - current->start_index+1));
                            int origin_left_delta_bit = current->estimate_bit;
                            int origin_right_delta_bit = next->estimate_bit;
                            int origin_cost = (right_index - next->start_index + 1) * origin_right_delta_bit + (next->start_index - current->start_index) * origin_left_delta_bit;
                            int merged_cost = new_bit * (right_index - current->start_index + 1);

                            if (merged_cost - origin_cost < overhead)
                            {
                                // merge
                                current->end_index = right_index;
                                current->next = next->next;
                                next->next->prev = current;
                                current->double_delta_next = next->double_delta_next;
                                current->estimate_bit = new_bit;
                                delete next;
                                next = current->next;
                            }
                            else
                            {
                                break;
                            }
                        }

                        current = current->next;
                    }
                    else
                    {
                        current = current->next;
                    }
                }
            }

            double end_timer = getNow();
            split_time += (end_timer - start_timer);

            current = (&head)->next;
            while (current->next != &tail)
            {
                segment_index.push_back(current->start_index);
                // std::cout << current->start_index << " " << current->end_index << " " << current->estimate_bit << std::endl;
                current = current->next;
            }
            if(nvalue==70){
                std::cout<<segment_index.size()<<std::endl;
            }
            int segment_total = segment_index.size();
            segment_index.push_back(length);

            total_byte = 0;
            for (int i = 0; i < segment_total; i++)
            {
                uint64_t tmp_size = newsegment_size(segment_index[i], segment_index[i + 1] - 1);
                total_byte += tmp_size;
                // std::cout<<segment_index[i]<<" "<<segment_index[i + 1] - 1<<" "<<total_byte<<std::endl;
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
                merge(nvalue,length);
                // merge_both_direction(nvalue, length);
                double compressrate = (total_byte) * 100.0 / (sizeof(T) * block_size * 1.0);
                // std::cout << "try " << iter << " segment number " << (int)segment_index.size() << " resulting compression rate: " << std::setprecision(4) << compressrate << std::endl;
                cost_decline = cost_decline - total_byte;
                double cost_decline_percent = cost_decline * 100.0 / (sizeof(T) * block_size * 1.0);
                if (cost_decline_percent < 0.01) {
                    break;
                }

            }

            int segment_number = (int)segment_index.size();
            if (segment_number <= 1)
            {
                segment_index.clear();
                segment_index.push_back(0);
                segment_index.push_back(length);
                newsegment(0,length - 1, nvalue);
                // std::cout<<nvalue<<" "<<segment_index[1]-segment_index[0]<<" "<< start_key[0]<<" "<< slope[0]<<std::endl;
            }
            else
            {
                segment_index.push_back(length);
                for (int i = 0; i < segment_number; i++)
                {
                    newsegment(segment_index[i], segment_index[i + 1] - 1, nvalue);
                    // std::cout<<nvalue<< " "<<segment_index[i]<<" "<<segment_index[i+1]-1<<std::endl;
                    // std::cout<<nvalue<<" "<<segment_index[i + 1]-segment_index[i]<<" "<< start_key[i]<<" "<< slope[i]<<std::endl;
                }
            }
            segment_index.pop_back();

            end_timer = getNow();
            merge_time += (end_timer - start_timer);

            for (auto item : segment_index)
            {
                art_build_vec.push_back((KeyValue<uint32_t>){nvalue*block_size + item, segment_index_total_idx});
                // std::cout<<item<<std::endl;
                // auto tmp = btree_total.insert(std::make_pair(item, segment_index_total_idx));
                // auto tmp = alex_tree.insert(item, segment_index_total_idx);
                alex_build_vec.push_back(std::make_pair(nvalue*block_size +item, segment_index_total_idx));
                segment_index_total_idx++;
            }

            if (nvalue == block_num - 1)
            {
                // auto tmp = btree_total.insert(std::make_pair(block_num * block_size, segment_index_total_idx));
                // std::cout<<block_num * block_size<<" "<<segment_index_total_idx<<std::endl;
                // auto tmp = alex_tree.insert(block_num * block_size, segment_index_total_idx);
                art_build_vec.push_back((KeyValue<uint32_t>){block_num * block_size, segment_index_total_idx});
                alex_build_vec.push_back(std::make_pair(block_num * block_size, segment_index_total_idx));
            }
            segment_index.clear();
            segment_length.clear();
            total_byte = 0;
            current = (&head)->next;
            Segment_pol<int64_t> *next = current->next;
            while (current != &tail && current->next != &tail)
            {
                next = current->next;
                delete current;
                current = next;
            }
            (&head)->next = &tail;
            (&tail)->prev = &head;
            // for(int i=0; i< segment_index_total.size()-1; i++){
                
            //     std::cout<<segment_index_total[i]<<" "<< (double)newsegment_size(segment_index_total[i], segment_index_total[i+1]-1) / (double)((segment_index_total[i+1] - segment_index_total[i])*8)<<std::endl;
            // }

            return res;
        }

        void merge(int nvalue, int length) {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = segment_index[0]; // before the start_index is the finished blocks
            int segment_num = 0; // the current segment index
            int newsegment_num = 0;
            int total_segments = segment_index.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;
            segment_index.push_back(length);
     
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



        void merge_both_direction(int nvalue, int length)
        {
            // this function is to merge blocks in block_start_vec to large blocks
            int start_index = segment_index[0];        // before the start_index is the finished blocks
            int segment_num = 0;                       // the current segment index
            int total_segments = segment_index.size(); // the total number of segments
            uint64_t totalbyte_after_merge = 0;

            segment_index.push_back(length);
            std::vector<uint32_t> new_segment_index;
            std::vector<uint32_t> new_segment_length;
            new_segment_index.emplace_back(segment_index[segment_num]);
            new_segment_length.emplace_back(segment_length[segment_num]);
            totalbyte_after_merge += segment_length[segment_num];
            segment_num++;

            while (segment_num < total_segments)
            {
                // std::cout<<"segment_num: "<<segment_num <<" / "<<total_segments<<std::endl;

                if (segment_num == total_segments - 1)
                {
                    // only can try merging with former one
                    int last_merged_segment = new_segment_length.size() - 1;
                    uint32_t init_cost_front = segment_length[segment_num] + new_segment_length[last_merged_segment];
                    uint32_t merge_cost_front = newsegment_size(new_segment_index[last_merged_segment], segment_index[segment_num + 1] - 1);
                    int saved_cost_front = init_cost_front - merge_cost_front;
                    if (saved_cost_front > 0)
                    {
                        totalbyte_after_merge -= new_segment_length[new_segment_length.size() - 1];
                        new_segment_length[new_segment_length.size() - 1] = merge_cost_front;
                        totalbyte_after_merge += merge_cost_front;
                        segment_num++;
                    }
                    else
                    {
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
                if (saved_cost <= 0)
                {
                    // do not merge
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(segment_length[segment_num]);
                    totalbyte_after_merge += segment_length[segment_num];
                    start_index = segment_index[segment_num + 1];
                    segment_num++;

                    continue;
                }
                if (saved_cost_back > saved_cost_front)
                {
                    // merge with back
                    new_segment_index.emplace_back(segment_index[segment_num]);
                    new_segment_length.emplace_back(merge_cost_back);
                    totalbyte_after_merge += merge_cost_back;
                    start_index = segment_index[segment_num + 2];
                    segment_num += 2;
                    // std::cout<<segment_num<<std::endl;
                }
                else
                {
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

        T *decodeArray8(const size_t length, T *out, size_t nvalue)
        {
            T *res = out;
            // start_index + bit + theta0 + theta1 + numbers + delta
            segment_index_total.push_back(length);
            std::array<double, degree + 1> theta;
            uint8_t maxerror;
            for (int i = 0; i < block_start_vec_total.size(); i++)
            {
                int segment_length = segment_index_total[i + 1] - segment_index_total[i];
                uint8_t *tmpin = block_start_vec_total[i];
                // debug
                int start_ind = 0;
                memcpy(&start_ind, tmpin, sizeof(int));

                // debug
                tmpin += sizeof(uint32_t);
                maxerror = tmpin[0];
                tmpin++;
                if (maxerror == 255)
                {
                    T tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }
                if (maxerror == 254)
                {
                    T tmp_val;
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    memcpy(&tmp_val, tmpin + sizeof(T), sizeof(tmp_val));
                    res[0] = tmp_val;
                    res++;
                    continue;
                }
                if (maxerror == sizeof(T) * 8)
                {
                    memcpy(res, tmpin, sizeof(T) * segment_length);
                    res += segment_length;
                    continue;
                }

                for (auto i = 0; i < degree + 1; i++)
                {
                    memcpy(&theta[i], tmpin, sizeof(double));
                    tmpin += sizeof(double);
                }
                andviane::Polynomial<degree> poly(theta);

                if (maxerror)
                {
                    read_all_bit_fix_poly<T,degree,double>(tmpin, 0, 0, segment_length, maxerror, res, &poly);
                }
                else
                {
                    for (int j = 0; j < segment_length; j++)
                    {
                        res[j] = (long long)(poly(j));
                    }
                }
                res += segment_length;
            }
            return out;
        }

        int get_segment_id(int to_find)
        {
            // int segment_id = art.upper_bound_new(to_find, search_node) - 1;
            // int segment_id = lower_bound(to_find, segment_index_total.size(), segment_index_total);
            int segment_id = alex_tree.upper_bound(to_find).payload() - 1;
            __builtin_prefetch(block_start_vec_total.data() + segment_id, 0, 3);
            return segment_id;
        }

        T randomdecodeArray8(int segment_id, uint8_t *in, int to_find, uint32_t *out, size_t nvalue)
        {

            uint8_t *this_block = block_start_vec_total[segment_id];

            uint8_t *tmpin = this_block;

            uint32_t start_ind;
            memcpy(&start_ind, tmpin, 4);
            tmpin += 4;

            uint8_t maxerror;
            memcpy(&maxerror, tmpin, 1);
            tmpin++;
            if (maxerror == sizeof(T) * 8)
            {
                T tmp_val = reinterpret_cast<T *>(tmpin)[to_find - start_ind];
                return tmp_val;
            }

            T tmp_val = 0;
            if (maxerror == 255)
            {
                memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                return tmp_val;
            }
            if (maxerror == 254)
            {
                if (to_find - start_ind == 0)
                {
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                }
                else
                {
                    tmpin += sizeof(tmp_val);
                    memcpy(&tmp_val, tmpin, sizeof(tmp_val));
                }
                return tmp_val;
            }


            std::array<double, degree + 1> theta;
            for (auto i = 0; i < degree + 1; i++)
            {
                memcpy(&theta[i], tmpin, sizeof(double));
                tmpin += sizeof(double);
            }
            andviane::Polynomial<degree> poly(theta);
            if (maxerror)
            {
                tmp_val = read_bit_fix_int_wo_round_poly<T,degree,double>(tmpin, maxerror, to_find-start_ind, &poly);
            }
            else
            {
                tmp_val = (int64_t)poly(to_find-start_ind);
            }

            return tmp_val;
        }

        uint64_t summation(uint8_t *in, const size_t l, size_t nvalue)
        {

            return 0;
        }
        uint32_t *encodeArray(uint32_t *in, const size_t length, uint32_t *out,
                              size_t nvalue)
        {
            std::cout << "Haven't implement. Please try uint8_t one..." << std::endl;
            return out;
        }
        uint32_t *decodeArray(uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue)
        {
            std::cout << "Haven't implement. Please try uint8_t one..." << std::endl;
            return out;
        }
        uint32_t randomdecodeArray(uint32_t *in, const size_t l, uint32_t *out, size_t nvalue)
        {
            std::cout << "Haven't implement. Please try uint8_t one..." << std::endl;
            return 1;
        }
        uint32_t get_total_byte()
        {
            return total_byte_total;
        }
        uint32_t get_total_blocks()
        {
            // std::cout << "split time " << split_time << std::endl;
            // std::cout << "merge time " << merge_time << std::endl;
            return block_start_vec_total.size();
        }
    };

} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
