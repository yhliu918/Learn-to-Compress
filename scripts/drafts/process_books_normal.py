import numpy as np
# filename = '/home/lyh/Learn-to-Compress/integer_data/books_200M_uint32.txt'
# ff = open(filename, 'r')
# data = [int(item) for item in ff.readlines()]

f = open("/home/lyh/Learn-to-Compress/integer_data/fb_200M_uint64", "r")
data = np.fromfile(f, dtype=np.uint64)[1:]
print(len(data))
start = 0
end = 744520
with open('/home/lyh/Learn-to-Compress/integer_data/fb_part.txt','w') as fout:
    for i in range(start,end):
        fout.write(str(data[i])+'\n')
