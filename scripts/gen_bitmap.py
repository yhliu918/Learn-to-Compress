import numpy as np
import struct
from scipy.stats import norm, lognorm
import os

# any arbitrary seed value will do, but this one is clearly the best.
np.random.seed(seed=42) 

NUM_KEYS = 289000


#prob = [0.00001, 0.0001, 0.0005, 0.001, 0.01]
prob = [0.1]
# 0.001%, 0.01%, 0.05%, 0.1%, 1%
for prob_ in prob:
    print("Generating random bitmap...")
    print(prob_)
    file_name = "../data/bitmap_random/bitmap_random_"+str(prob_)+"_"+str(NUM_KEYS)+".txt"
    if not os.path.exists(file_name):
        count = 0
        keys = np.random.choice([0, 1], size=NUM_KEYS, p=[1.0-prob_, prob_])
        for i in range(NUM_KEYS):
            if keys[i] == 1:
                count+=1
        print(count)
        keys = np.array(keys)
        keys = keys.astype(np.uint32)
    

        np.savetxt(file_name, keys,fmt='%d')



'''
selectivity_ = [ 0.001, 0.01] 
cluster_ = [10,100, 1000]
print("Generating cluster bitmap...")
for selectivity in selectivity_:
    for cluster in cluster_:
        if cluster > selectivity*NUM_KEYS:
            continue
        
        file_name = "../data/bitmap_cluster/bitmap_cluster_"+str(selectivity)+"_"+str(cluster)+"_"+str(NUM_KEYS)+".txt"
        print(file_name)
        if not os.path.exists(file_name):
            start_ind = []
            keys = []
            # selectivity * NUM_KEYS
            count = 0
            for i in range(cluster):
                ind = np.random.randint(0,(1-selectivity)*NUM_KEYS/cluster) #avoid that the last 1 is in the next seg
                start_ind.append(i*NUM_KEYS/cluster + ind)
            for i in range(cluster):
                for j in range(int(NUM_KEYS/cluster)):
                    tmp = i*NUM_KEYS/cluster + j
                    if tmp>=start_ind[i] and tmp<start_ind[i]+selectivity*NUM_KEYS/cluster:
                        keys.append(1)
                        count+=1
                    else:
                        keys.append(0)
            print(count)
            keys = np.array(keys)
            keys = keys.astype(np.uint32)
    

            np.savetxt(file_name, keys,fmt='%d')
'''