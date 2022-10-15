
import random

N = 10000
for iter in range(1000):
    d = {}
    a = random.uniform(0, 1)
    b = random.uniform(0, 1)
    print(a, b)
    pred = b
    groundtruth = 0
    for i in range(N):
        groundtruth+=int(b + a*i)

    our_result = 0
    origin_k = int((1-b)/a)
    for i in range(origin_k):
        d.setdefault(0,[]).append(i)

    value = 1
    step = int(1/a)
    epsilon = 1 - step*a
    start_idx = origin_k+1
    end_idx = 0

    residual = b + a + origin_k*a - 1

    pred = b
    i = start_idx
    while i<N:
        end_idx = start_idx + step
        
        if residual - epsilon >=0:
            for j in range(start_idx, end_idx):
                d.setdefault(value,[]).append(j)
                i+=1
                if(i>N-1):
                    break
            residual -= epsilon
            value += 1
            start_idx = end_idx
        else:
            for j in range(start_idx, end_idx+1):
                d.setdefault(value,[]).append(j)
                i+=1
                if(i>N-1):
                    break
            residual = residual - epsilon + a
            value += 1
            start_idx = end_idx+1

        # print(start_idx, end_idx,i)
    # print(d)
    for key in d.keys():
        our_result += int(key) * len(d[key])

    print(groundtruth, our_result)
    if(groundtruth != our_result):
        print("error")
        break

        

    