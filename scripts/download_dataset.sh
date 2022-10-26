#linear, normal generated from gen_norm.py
python3 gen_norm.py
#books
wget -O books https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/5YTV8K
zstd -d books -o books_200M_uint32
#poisson generated from poisson_randomdie.py
python3 poisson_randomdie.py
#ml
wget -O https://archive.ics.uci.edu/ml/machine-learning-databases/00515/data.zip
#fb, wiki download from sosd
wget -O fb  https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/EATHF7
zstd -d fb -o fb_200M_uint64
wget -O wiki https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/JGVF9A/SVN8PI
zstd -d wiki -o wiki_200M_uint64
