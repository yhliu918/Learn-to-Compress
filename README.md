# Learn-to-Compress
To install dependent packages, run the following script if you have root access:
```
./scripts/setup_dependencies.sh
```

To have a quick start, you can run following commands to download and generate data needed.
```
mkdir data && cd data
wget -O ori_file https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/5YTV8K
zstd -d ori_file -o books_200M_uint32
rm ori_file
cd ../
python ./scripts/gen_norm.py
```
Then, you can run the command (Need -std=c++11 supported)
```
mkdir build && cd build
cmake ..
make
./example
```
You can change the API of different compression methods in the example.cpp
