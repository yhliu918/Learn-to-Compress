import os
import sys

datasets = ["linear_200M_uint32.txt","normal_200M_uint32.txt","books_200M_uint32.txt","poisson_randomdie.txt", "ml_timestamp.txt","house_price.txt", "movieid.txt","fb_200M_uint64","wiki_200M_uint64"]
# datasets = [ "fb/fb-289000.txt", "wf/wiki.txt", "wf/newman.txt", "unigram_freq_count.txt", "hu_freq.txt", "house_price.txt", "movieid.txt"]
method = ["delta_int","leco_int"]
#"leco_int_double","leco_int_double_wo_round"
# method = ["leco_int"]
model_size = {
    "delta_int":4,
    "leco_int":8,
    "leco_int_double":12,
    "leco_int_double_wo_round":12
}
blocks = {
    "linear_200M_uint32.txt": 100000,
    "normal_200M_uint32.txt": 100000,
    "books_200M_uint32.txt":100000, 
    "poisson_randomdie.txt":100000,
    "ml_timestamp.txt":10000,
    "house_price.txt":1, 
    "movieid.txt":1000,
    "fb_200M_uint64":100000,
    "wiki_200M_uint64":100000
}

delta = {
    "linear_200M_uint32.txt": 2,
    "normal_200M_uint32.txt": 2,
    "books_200M_uint32.txt":6, 
    "poisson_randomdie.txt":12,
    "ml_timestamp.txt":6,
    "house_price.txt":2, 
    "movieid.txt":12,
    "fb_200M_uint64":12,
    "wiki_200M_uint64":0
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

for dataset in datasets:
    for codec in method:
        if dataset == 'linear_200M_uint32.txt' or dataset == 'normal_200M_uint32.txt':
            if codec == 'leco_int':
                codec = 'leco_int_double'
                
        if dataset in ['poisson_randomdie.txt', 'ml_timestamp.txt','fb_200M_uint64','wiki_200M_uint64']:
            if codec == 'leco_int':
                codec = 'leco_int_double_wo_round'
        delta_use = delta[dataset]
        if codec =='delta_int' :
            delta_use = 4
        print(" ./{} {} {} {} {} {}".format(codec, dataset, blocks[dataset],delta_use, model_size[codec], binary[dataset]))