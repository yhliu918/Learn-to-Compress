import os
import sys

datasets = ["linear_200M_uint32.txt","normal_200M_uint32.txt","books_200M_uint32.txt","house_price.txt", "movieid.txt","poisson_randomdie.txt", "ml_timestamp.txt","fb_200M_uint64","wiki_200M_uint64"]
method = ["FOR_my","delta_my", "piecewise_fix_op_max"]
# method = [ "piecewise_fix_op","piecewise_fix_op_max"]
model_size = {
    "FOR_my":4,
    "delta_my":4,
    "piecewise_fix_op":16,
    "piecewise_fix_op_max":12,
    "FOR_int_fix":4,
    "Delta_int_fix":4,
    "leco_int_fix":12,
}
blocks = {
    "linear_200M_uint32.txt": 100000,
    "normal_200M_uint32.txt": 100000,
    "books_200M_uint32.txt":1000000, 
    "house_price.txt":1000, 
    "movieid.txt":100000,
    "poisson_randomdie.txt":1000000,
    "ml_timestamp.txt":100000,
    "fb_200M_uint64":1000000,
    "wiki_200M_uint64":1000000
}
binary ={
    "linear_200M_uint32.txt": 0,
    "normal_200M_uint32.txt": 0,
    "books_200M_uint32.txt":0, 
    "poisson_randomdie.txt":0,
    "ml_timestamp.txt":0,
    "house_price.txt":0, 
    "movieid.txt":0,
    "fb_200M_uint64":1,
    "wiki_200M_uint64":1
}
big_int ={
    "linear_200M_uint32.txt": 0,
    "normal_200M_uint32.txt": 0,
    "books_200M_uint32.txt":0, 
    "poisson_randomdie.txt":1,
    "ml_timestamp.txt":1,
    "house_price.txt":0, 
    "movieid.txt":0,
    "fb_200M_uint64":1,
    "wiki_200M_uint64":1
}

for dataset in datasets:
    if big_int[dataset]:
        for codec in method:
            if codec == 'FOR_my':
                codec = 'FOR_int_fix'
            elif codec == 'delta_my':
                codec = 'Delta_int_fix'
            elif codec == 'piecewise_fix_op_max':
                codec = 'leco_int_fix'
            print(" ./{} {} {} {} {}".format(codec, dataset, blocks[dataset], model_size[codec], binary[dataset]))
    else:
        for codec in method:
            print(" ./fix_int {} {} {} {}".format(codec, dataset, blocks[dataset], model_size[codec]))