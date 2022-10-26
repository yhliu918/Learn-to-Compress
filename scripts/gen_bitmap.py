import numpy as np
import struct
from scipy.stats import norm, lognorm
import os
import threading

# any arbitrary seed value will do, but this one is clearly the best.
np.random.seed(seed=42) 

# NUM_KEYS = 200000000
# NUM_KEYS = 289000
# NUM_KEYS_list = [1000000]
# NUM_KEYS_list = [200000000, 199999994]
# NUM_KEYS_list = [14057565]
NUM_KEYS_list = [200015910]
# NUM_KEYS_list = [200000000, 199999994, 289000, 2076000, 233000]


prob = [0.00001, 0.0001, 0.0005, 0.001, 0.01, 0.1]
# prob = [0.1]
# 0.001%, 0.01%, 0.05%, 0.1%, 1%

def uniform():
    for NUM_KEYS in NUM_KEYS_list:
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

def myThread(divide_base, num_cluster, prob_, NUM_KEYS):
    factor = 1 / divide_base**num_cluster
    print("Generating random bitmap...")
    print(prob_)
    file_name = "../data/bitmap_random_cluster/bitmap_random_"+str(prob_)+"_"+str(NUM_KEYS)+".txt"
    # if not os.path.exists(file_name):
    count = 0
    keys_list = []
    sum_prob = 0
    for i in range(num_cluster - 1):
        keys_list.append(np.random.choice([0, 1], size=NUM_KEYS//num_cluster, p=[1.0-prob_*factor, prob_*factor]))
        sum_prob += factor
        factor *= divide_base
    keys_list.append(np.random.choice([0, 1], size=NUM_KEYS//num_cluster, p=[1.0-prob_*(num_cluster - sum_prob), prob_*(num_cluster - sum_prob)]))
    for i in range(NUM_KEYS//num_cluster):
        for j in keys_list:
            if j[i] == 1:
                count += 1
    print(count)
    keys = np.concatenate(tuple(keys_list))
    keys = keys.astype(np.uint32)
    np.savetxt(file_name, keys,fmt='%d')
    
def cluster():
    num_cluster = 10
    divide_base = 2
    ts = []
    for NUM_KEYS in NUM_KEYS_list:
        for prob_ in prob:
            t = threading.Thread(target=myThread, args=(divide_base, num_cluster, prob_, NUM_KEYS))
            t.start()
            ts.append(t)
    for t in ts:
        t.join()


def cluster_v0():
    for NUM_KEYS in NUM_KEYS_list:
        selectivity_ = prob
        cluster_ = [10]
        print("Generating cluster bitmap...")
        for selectivity in selectivity_:
            for cluster in cluster_:
                if cluster > selectivity*NUM_KEYS:
                    continue
                
                file_name = f"../data/bitmap/bitmap_random_{selectivity}_{NUM_KEYS}.txt"
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

if __name__ == "__main__":
    cluster()
    # cluster_v0()