import numpy as np
import struct
from scipy.stats import norm, lognorm
import os

# any arbitrary seed value will do, but this one is clearly the best.
np.random.seed(seed=42) 

NUM_KEYS = 200000000

# print("Generating pieceiwse linear data...")
# if not os.path.exists("../data/standard/piecewise_linear_200M_uint32.txt"):
#     print("32 bit...")
#     count = 0
#     keys=[]
#     last = 0
#     while count<NUM_KEYS:
#         length = np.random.randint(0,1000)
#         slope = np.random.random_sample()+np.random.randint(30)
#         print(length,slope)
#         for i in range(length):
#             tmp = last+slope
#             keys.append(tmp)
#             last = tmp
#             count+=1
#             if count>=NUM_KEYS:
#                 break
#     keys = np.array(keys)
#     keys = keys.astype(np.uint32)
    

#     np.savetxt("../data/standard/piecewise_linear_200M_uint32.txt", keys,fmt='%d')

print("Generating linear data...")
if not os.path.exists("data/linear_200M_uint32.txt"):
    print("32 bit...")
    keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]
    keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
    keys *= 2**32 - 1
    keys = keys.astype(np.uint32)
    np.savetxt("data/linear_200M_uint32.txt", keys,fmt='%d')

'''
print("Generating noisy linear data...")
if not os.path.exists("../data/test.txt"):
    print("32 bit...")
    keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]
    noise=np.random.normal(0, 1000000, int(NUM_KEYS/100) )
    noise[0]=0
    noise[-1]=0
    k=0
    keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
    keys *= 2**32 - 1
    for i in range(NUM_KEYS):
        if np.random.randint(0,999)==0:
            keys[i]=keys[i]+noise[k]
            if keys[i]<0:
                keys[i]=0
            elif keys[i]>2**32-1:
                keys[i]=2**32-1
            print(str(i)+' '+str(k)+' '+str(keys[i]))
            k+=1
    keys = keys.astype(np.uint32)
    

    np.savetxt("../data/test.txt", keys,fmt='%d')
'''

# print("Generating noisy normal data...")   
# if not os.path.exists("data/noisynormal_200M_uint32.txt"):
#     print("32 bit...")
#     keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]
#     noise=np.random.normal(0, 5, NUM_KEYS )
#     noise[0]=0
#     noise[-1]=0
#     # for some reason, the PPF function seems to use quadratic memory
#     # with the size of its input.
#     keys = np.array_split(keys, 1000)
#     keys = [norm.ppf(x) for x in keys]
#     keys = np.array(keys).flatten()
    
#     keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
#     keys *= 2**32 - 1
#     keys+=noise
#     keys = keys.astype(np.uint32)
#     np.savetxt("data/noisynormal_200M_uint32.txt", keys,fmt='%d')

print("Generating normal data...")

if not os.path.exists("data/normal_200M_uint32.txt"):
    print("32 bit...")
    keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]

    # for some reason, the PPF function seems to use quadratic memory
    # with the size of its input.
    keys = np.array_split(keys, 1000)
    keys = [norm.ppf(x) for x in keys]
    keys = np.array(keys).flatten()

    keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
    keys *= 2**32 - 1
    keys = keys.astype(np.uint32)
    np.savetxt("data/normal_200M_uint32.txt", keys,fmt='%d')

# print("Generating noisy normal data...")   
# if not os.path.exists("data/noisynormal_200M_uint32.txt"):
#     print("32 bit...")
#     keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]
#     noise=np.random.normal(0, 5, NUM_KEYS )
#     noise[0]=0
#     noise[-1]=0
#     # for some reason, the PPF function seems to use quadratic memory
#     # with the size of its input.
#     keys = np.array_split(keys, 1000)
#     keys = [norm.ppf(x) for x in keys]
#     keys = np.array(keys).flatten()
    
#     keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
#     keys *= 2**32 - 1
#     keys+=noise
#     keys = keys.astype(np.uint32)
#     np.savetxt("data/noisynormal_200M_uint32.txt", keys,fmt='%d')



# print("Generating log normal data...")
# if not os.path.exists("data/lognormal_200M_uint32.txt"):
#     print("32 bit...")
#     keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]

#     # using a sigma of 2 for the 32 bit keys produces WAY too many
#     # duplicates, so we will deviate from the RMI paper
#     # and use 1.
#     keys = np.array_split(keys, 1000)
#     keys = [lognorm.ppf(x, 1) for x in keys]
#     keys = np.array(keys).flatten()

#     keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
#     keys *= 2**32 - 1
#     keys = keys.astype(np.uint32)
#     np.savetxt("data/lognormal_200M_uint32.txt", keys,fmt='%d')


# print("Generating books data...")
# if not os.path.exists("data/books_200M_uint32.txt"):
#     print("32 bit...")
#     import chardet
#     keys = np.fromfile("data/books_200M_uint32", dtype=np.uint32)[2:]
#     np.savetxt("data/books_200M_uint32.txt", keys,fmt='%d')   
