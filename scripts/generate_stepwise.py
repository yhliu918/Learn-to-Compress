import random
interval = 5
block_size = 10000
N = 100000000
res = 3
l = []
for j in range(int(N/block_size)):
    for i in range(block_size):
        l.append(res+interval*i+random.randint(-3, 3))

with open('noisy_stepwise.txt','w')as ff:
    for item in l:
        ff.write(str(item)+'\n')