import sys, os
datasets = ["fb/fb-289000.txt","wf/wiki.txt","wf/newman.txt", "house_price.txt","movieid.txt"]
logname = ['fb', 'wiki','newman','house_price','movieid']
for i in range(len(datasets)):
    print("python3 leco_lp_cost.py {} {} > {}.log  ".format(datasets[i],logname[i], logname[i]))