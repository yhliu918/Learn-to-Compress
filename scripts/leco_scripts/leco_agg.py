
import random

N = 200
for iter in range(1):

    # a = random.uniform(0, 1)
    # b = random.uniform(0, 1)
    # print(a, b)
    a = 0.037505500999998276
    b = -0.9017909989997861
    print(a, b)
    pred = b
    groundtruth = 0
    for i in range(N):
        groundtruth+=int(b + a*i)

    our_result = 0
    value_number_list = []
    origin_k = int((1-b)/a)
    value_number_list.append(origin_k)
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
            number = end_idx - start_idx
            i += number 
            if i>N-1:
                value_number_list.append(N-start_idx)
                break
            value_number_list.append(number)
            residual -= epsilon
            value += 1
            start_idx = end_idx
        else:
            number = end_idx - start_idx +1
            i += number 
            if i>N-1:
                value_number_list.append(N-start_idx)
                break
            value_number_list.append(number)

            residual = residual - epsilon + a
            value += 1
            start_idx = end_idx+1


    for key in range(len(value_number_list)):
        our_result += value_number_list[key] * key
    print(groundtruth, our_result)
    if(groundtruth != our_result):
        print("error")
        break

        

    