#ifndef RANK_H_
#define RANK_H_


#include "popcount.h"
/*
class Rank{
public:
    int* rank_lut_;
    int basic_block_size_=64;
void set_basic_block_size(int basic_block){
    basic_block_size_ = basic_block;
}
int rank( int pos, uint64_t *bits_  ){//basic_block_size_ is like 64

    int block_id = pos/basic_block_size_;
    int offset = pos & (basic_block_size_ - 1);
    //num_bits_ is the block size
    return (rank_lut_[block_id] + popcountLinear(bits_, block_id , offset + 1));
}

void initRankLut(int num_bits_, uint64_t *bits_  ) {
        int num_blocks = num_bits_ / basic_block_size_ + 1;
        rank_lut_ = new int[num_blocks];

        int cumu_rank = 0;
        for (int i = 0; i < num_blocks - 1; i++) {
            rank_lut_[i] = cumu_rank;
            cumu_rank += popcountLinear(bits_, i, basic_block_size_);
        }
        rank_lut_[num_blocks - 1] = cumu_rank;
}

};
*/

int rank( int pos, uint64_t *bits_ ,int*rank_lut_,int basic_block_size_ =64){//basic_block_size_ is like 64

    int block_id = pos/basic_block_size_;
    int offset = pos & (basic_block_size_ - 1);
    //num_bits_ is the block size
    return (rank_lut_[block_id] + popcountLinear(bits_, block_id , offset + 1));
}

void initRankLut(int num_bits_, uint64_t *bits_ ,int*rank_lut_ ,int basic_block_size_=64) {
        int num_blocks = ceil((double)num_bits_ / (double)basic_block_size_ );
        int cumu_rank = 0;
        for (int i = 0; i < num_blocks - 1; i++) {
            rank_lut_[i] = cumu_rank;
            cumu_rank += popcountLinear(bits_, i, basic_block_size_);
        }
        rank_lut_[num_blocks - 1] = cumu_rank;

}


#endif // RANK_H_
