import os
import sys

datasets = ["books_200M_uint32.txt", "fb/fb-289000.txt", "wf/wiki.txt", "wf/newman.txt", "unigram_freq_count.txt", "hu_freq.txt", "house_price.txt", "movieid.txt"]

method = ["delta_int","leco_int"]
model_size = {
    "delta_int":4,
    "leco_int":8
}
blocks = {
    "books_200M_uint32.txt":100000, 
    "fb/fb-289000.txt":10, 
    "wf/wiki.txt":10,
    "wf/newman.txt":1,
    "unigram_freq_count.txt":10, 
    "hu_freq.txt":10, 
    "house_price.txt":1, 
    "movieid.txt":1000
}

delta = {
    "books_200M_uint32.txt":6, 
    "fb/fb-289000.txt":2, 
    "wf/wiki.txt":4,
    "wf/newman.txt":2,
    "unigram_freq_count.txt":0, 
    "hu_freq.txt":2, 
    "house_price.txt":2, 
    "movieid.txt":12
}

for dataset in datasets:
    for codec in method:
        os.system(" ./{} {} {} {} {}".format(codec, dataset, blocks[dataset],delta[dataset], model_size[codec]))