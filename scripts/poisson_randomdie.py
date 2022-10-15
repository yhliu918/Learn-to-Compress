import sys
import numpy as np
from random import randint
OUTER_GAP= 1000000000000
GAP = 2000000
#     200000000
#     2000000
RAND_SENSOR_DIE=500
NUM_SENSOR = 2000
# parquet param, commented off
# NUM_EVENT_PER_SENSOR = 114296
# OUTFILE_PATH = "../data/poisson_timestamps_EVENT_50000_SENSOR_2000_randomdie_OUTER_1000s_INNER_2ms_200M.csv"
NUM_EVENT_PER_SENSOR = 50000

OUTFILE_PATH = "../integer_data/poisson_timestamps_EVENT_50000_SENSOR_2000_randomdie_OUTER_1000s_INNER_2ms_100M.csv"

f_out = open(OUTFILE_PATH, 'w')

outer_gap_list = np.random.poisson(OUTER_GAP, NUM_EVENT_PER_SENSOR)
for j in range(1, NUM_EVENT_PER_SENSOR) :
    outer_gap_list[j] += outer_gap_list[j-1]
print(outer_gap_list)

for i in range(0, NUM_EVENT_PER_SENSOR) :
    random_sensor_num = randint(0, RAND_SENSOR_DIE)
    remain_sensor_num = NUM_SENSOR - random_sensor_num
    data = np.random.poisson(GAP, remain_sensor_num)
    data[0] = outer_gap_list[i]
    # print(data[0])
    for j in range(1, remain_sensor_num) :
        data[j] += data[j-1]

    for k in range(0, remain_sensor_num) :
        f_out.write(str(data[k]) + "\n")
        
f_out.close()