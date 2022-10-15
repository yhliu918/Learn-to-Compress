import sys
import numpy as np
from random import randint
OUTER_GAP= 100000000
GAP = 20000
#     200000000
#     2000000
NUM_SENSOR = 20
NUM_EVENT_PER_SENSOR = 500000

OUTFILE_PATH = "/home/lyh/Learn-to-Compress/integer_data/poisson_sparse.txt"

f_out = open(OUTFILE_PATH, 'w')

outer_gap_list = np.random.poisson(OUTER_GAP, NUM_EVENT_PER_SENSOR)
for j in range(1, NUM_EVENT_PER_SENSOR) :
    outer_gap_list[j] += outer_gap_list[j-1]
print(outer_gap_list)

for i in range(0, NUM_EVENT_PER_SENSOR) :
    data = np.random.poisson(GAP, NUM_SENSOR)
    data[0] = outer_gap_list[i]
    # print(data[0])
    for j in range(1, NUM_SENSOR) :
        data[j] += data[j-1]

    for k in range(0, NUM_SENSOR) :
        f_out.write(str(data[k]) + "\n")
        
f_out.close()