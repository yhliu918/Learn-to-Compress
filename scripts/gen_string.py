import numpy as np
import struct
from scipy.stats import norm, lognorm
import os
import random
import pickle
np.random.seed(seed=42) 

NUM_KEYS = 1000000

def converttoascii(x):
    str = ''
    while x>1:
        if x%128 == 10:
            str = chr(x%128+1)+str
        else:
            str = chr(x%128)+str
        x = x//256
    return str



def converttoint(x):
    code = 0
    while len(x)>0:
        code = code*(1<<8)+ord(x[0])
        x = x[1:]
    return code
STRING_MIN = 'aaaaaaaa'
STRING_MAX = 'zzzzzzzz'

def generate_rand_string():
    return ''.join(random.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(8))

# outfile = "../data/string_linear_1M.txt"
'''
with open(outfile, 'w') as f:
    l = []
    for i in range(NUM_KEYS):
        l.append(generate_rand_string())
    l = sorted(l)
    for item in range(len(l)):
        # print(l[item])
        f.write(l[item]+'\n')
'''        



# print("Generating linear data...")
# outfile = "../string_data/test_hex.txt"
# with open(outfile, 'w') as f:
#     keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]
#     keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
#     keys *= (2**63)
#     # keys += converttoint(STRING_MIN)
#     for item in range(len(keys)):
#         tmpstr = hex(int(keys[item]))
#         f.write(hex(int(keys[item]))[2:]+'\n')



# outfile = "../string_data/test.txt"
# with open(outfile, 'w') as f:
#     for i in range(26):
#         cha = chr(ord('a')+i)
#         f.write(cha*8+'\n')
    
print("Generating linear data...")
outfile = "../string_data/string_linear_1M.txt"
with open(outfile, 'w') as f:
    keys = np.linspace(0, 1, NUM_KEYS + 2)[1:-1]
    keys = (keys - np.min(keys)) / (np.max(keys) - np.min(keys))
    keys *= (converttoint(STRING_MAX) - converttoint(STRING_MIN))
    keys += converttoint(STRING_MIN)
    # f.write(int64(NUM_KEYS))
    # pickle.dump(NUM_KEYS, f)
    for item in range(len(keys)):
        # print(converttoascii(int(keys[item])))
        # pickle.dump(converttoascii(int(keys[item])), f)
        f.write(converttoascii(int(keys[item]))+'\n')