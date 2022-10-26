import scipy
import numpy
import pandas as pd
import numpy as np
# np.set_printoptions(suppress=True)
import math
from scipy.optimize import linprog
import warnings
warnings.filterwarnings('ignore')
import sys, os



def calculate_lp(filename, blocks):
    data = []
    with open('../integer_data/'+filename,'r') as ff:
        for lines in ff:
            if len(lines) > 1:
                data.append(int(lines[:-1]))
    
    N = len(data)

    block_size = int(N/blocks)
    blocks = int(N/block_size)
    if blocks * block_size < N:
        blocks += 1
    print(blocks)
    total_bytes = 0
    for i in range(blocks):
        start_ind = block_size*i
        end_ind = min(block_size*(i+1),N)
        block_length = end_ind - start_ind 
        para = []
        b = []

        c = [0,0,1]
        for j in range(start_ind,end_ind):
            para.append([j-start_ind,1,-1])
            para.append([start_ind-j,-1,-1])
            b.append(data[j])
            b.append(-data[j])
        a_bound = (None, None)
        b_bound = (None, None)
        f_bound = (0, None)
        # print(len(c),len(para),len(b))
        res = linprog(c, A_ub=para, b_ub=b, options={'maxiter': 40},bounds=[a_bound, b_bound,f_bound])
        print('{:.5f} {:.5f}'.format(res.x[1],res.x[0]))
        # theta0=res.x[1]
        # theta1=res.x[0]
        # maxerror = 0
        # for j in range(start_ind,end_ind):
        #     tmp = data[j]-(theta0+theta1*j)
        #     if tmp > maxerror:
        #         maxerror = tmp
        # print('maxerror {}'.format(maxerror))
        # if i % 10 == 0:
        #     print(i)
        if res.x[-1] >= 0.1:
            res.x[-1] = math.ceil(res.x[-1])
            total_bytes+= math.ceil((math.log2(res.x[-1]+1)+1)*block_length /8)
    print('*'*20)
    print((total_bytes+blocks*9)/(N*4))
    return total_bytes+blocks*9
            


blocks = {
    "fb/fb-289000.txt":1000, 
    "wf/wiki.txt":10000,
    "wf/newman.txt":1000,
    "house_price.txt":1000, 
    "movieid.txt":100000,
    "books_200M_uint32.txt":1000000, 
    "fb_part.txt":37227,
}


# wf = open('./leco-lp.log','a')
# for name in datasets:
#     wf.write(name+'\n')
#     wf.write(str(calculate_lp(name,blocks[name]))+'\n')
name = sys.argv[1]
calculate_lp(name,blocks[name])