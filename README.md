# Learn-to-Compress
To have a quick start, you can run to download and generate data needed
```
wget -O - https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/5YTV8K | zstd -d > data/books_uint32_200M
python gen_norm.py

```
Then, you can run the command
```
g++ piecewise.cpp -o piecewise
./piecewise -f filename -e max_error
```
In which there are two hyper-parameter, filename in {linear,noisylinear,normal,noisynormal,lognormal,books}
and max_error is an arbitrary int you choose (2^k - 1 like 7,15,31,63,127... is recommended)
