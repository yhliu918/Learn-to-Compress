import sys, os
datasets = ["wf/wiki.txt", "house_price.txt","movieid.txt","books_200M_uint32.txt"]
logname = ['wiki','house_price','movieid','books_200M_uint32']
for i in range(len(datasets)):
    os.system("python3 leco_lp.py {} > {}.log  ".format(datasets[i],logname[i]))