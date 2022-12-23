import os
import sys


method = ["FOR_my","delta_my", "piecewise_fix_op_max_test"]
model_size = {
    "FOR_my":4,
    "delta_my":4,
    "piecewise_fix_op_max_test":12
}
datasets= ['/root/Learn-to-Compress/data/noisy_stepwise_10000.txt']
filter_val = {
    '/root/Learn-to-Compress/data/noisy_stepwise_10000.txt':[10000, 20000, 45000, 48000, 49900, 50000]
}
block={
    '/root/Learn-to-Compress/data/noisy_stepwise_10000.txt':100000,
}

for dataset in datasets:
    for bar in filter_val[dataset]:
        for codec in method:
            if codec == 'piecewise_fix_op_max_test':
                os.system(" ./leco_filter {} {} {} {}".format(dataset, block[dataset], model_size[codec], bar))
            else:
                os.system(" ./filter_test {} {} {} {} {}".format(codec, dataset,  block[dataset], model_size[codec], bar))