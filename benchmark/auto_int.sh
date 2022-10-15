#! /bin/bash
./delta_int linear_200M_uint32.txt 100000 4 4 0
./leco_int_double linear_200M_uint32.txt 100000 2 12 0
./delta_int normal_200M_uint32.txt 100000 4 4 0
./leco_int_double normal_200M_uint32.txt 100000 2 12 0
./delta_int books_200M_uint32.txt 100000 4 4 0
./leco_int books_200M_uint32.txt 100000 6 8 0
./delta_int poisson_randomdie.txt 100000 4 4 0
./leco_int_double_wo_round poisson_randomdie.txt 100000 12 12 0
./delta_int ml_timestamp.txt 10000 4 4 0
./leco_int_double_wo_round ml_timestamp.txt 10000 6 12 0
./delta_int house_price.txt 1 4 4 0
./leco_int house_price.txt 1 2 8 0
./delta_int movieid.txt 1000 4 4 0
./leco_int movieid.txt 1000 12 8 0
./delta_int fb_200M_uint64 100000 4 4 1
./leco_int_double_wo_round fb_200M_uint64 100000 12 12 1
./delta_int wiki_200M_uint64 100000 4 4 1
./leco_int_double_wo_round wiki_200M_uint64 100000 0 12 1