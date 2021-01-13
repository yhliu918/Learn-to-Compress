# Learn-to-Compress
To have a quick start, you can run to download and generate data needed
```
mkdir data

wget -O - https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/5YTV8K | zstd -d > data/books_uint32_200M

python gen_norm.py

```
Then, you can run the command (Need -std=c++11 supported)
```
mkdir build && cd build
cmake ..
make
```
You can change the API of different compression methods in the example.cpp
