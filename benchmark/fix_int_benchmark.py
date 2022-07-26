import os
import sys

datasets = ["linear_200M_uint32.txt", "normal_200M_uint32.txt","books_200M_uint32.txt", "fb/fb-289000.txt", "wf/wiki.txt", "wf/newman.txt", "unigram_freq_count.txt", "hu_freq.txt", "house_price.txt", "movieid.txt"]
method = ["delta_my","piecewise_fix_op_round", "piecewise_fix_op_max"]
model_size = {
    "delta_my":4,
    "piecewise_fix_op_round":16,
    "piecewise_fix_op_max":12
}
blocks = {
    "linear_200M_uint32.txt": 100000,
    "normal_200M_uint32.txt": 1000000,
    "books_200M_uint32.txt":1000000, 
    "fb/fb-289000.txt":1000, 
    "wf/wiki.txt":1000,
    "wf/newman.txt":1000,
    "unigram_freq_count.txt":1000, 
    "hu_freq.txt":1000, 
    "house_price.txt":1000, 
    "movieid.txt":100000
}

for dataset in datasets:
    for codec in method:
        os.system(" ./fix_int {} {} {} {}".format(codec, dataset, blocks[dataset], model_size[codec]))