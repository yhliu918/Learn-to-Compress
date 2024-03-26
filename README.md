# Leco: Learn-to-Compress
Related Research Paper:

**LeCo: Lightweight Compression via Learning Serial Correlations**

Author: Yihao Liu, Xinyu Zeng, Huanchen Zhang

[About to show up at **SIGMOD2024**]

```
@article{liu2023leco,
  title={LeCo: Lightweight Compression via Learning Serial Correlations},
  author={Liu, Yihao and Zeng, Xinyu and Zhang, Huanchen},
  journal={arXiv preprint arXiv:2306.15374},
  year={2023}
}
```

### Requirements:
Ubuntu 20.04.4 LTS

To install dependencies, run the following script:
```
sudo bash ./scripts/setup.sh
```

### Setup
To have a quick start, you can run following commands.
```
mkdir build && cd build
cmake ..
make -j16
```

Data sets used in the Evaluation part can be downloaded using `./script/download_dataset.sh`.

`house_price` and `movie_id` need to be manually downloaded from the following urls:
```
https://www.kaggle.com/datasets/ahmedshahriarsakib/usa-real-estate-dataset?select=realtor-dataset-100k.csv
https://www.kaggle.com/datasets/grouplens/movielens-20m-dataset?select=rating.csv
```
You can download and generate the data sets and store them in a directory `data` under the project root directory.

### Microbenchmark
The Microbenchmark of LeCo is located in `benchmark`, which can be run by:
```
python3 fix_int_benchmark.py > fix_int.sh
bash fix_int.sh > fix_int_benchmark_intel.log
python3 auto_int_benchmark.py > auto_int.sh
bash auto_int.sh > auto_int_benchmark_intel.log
```
To visualize the results, we provide the script in `scripts/paper_scripts/plot_cr_ra.ipynb`.
Except for system evaluation Figures, all other figures used in our paper is included in the above script.

### Layout
The file layout and explainations of our repository is as follows:
```
Learn-to-Compress
    \____ benchmark 
           \____ fix_int (benchmarking FOR, Delta_fix, LeCo_fix)
           \____ auto_int (benchmarking Delta_var, LeCo_var)
    \____ experiments (including all experiment of microbenchmark, as well as String Extensions)
    \____ headers (implementation of different compression schemes)
    \____ scripts (Setup/Data generation/Figure plot)
    \____ thirdparty
           \____ succinct (Elias-Fano)
           \____ FSST (String Baseline)
       
```
