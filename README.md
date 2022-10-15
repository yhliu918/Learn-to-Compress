# Leco: Learn-to-Compress

### Requirements:
```
Boost
gmp
Eigen3
```
To install dependent packages, run the following script if you have root access:
```
bash ./scripts/setup.sh
```

### Setup
To have a quick start, you can run following commands.
```
mkdir build && cd build
cmake ..
make -j16
```

All data sets used in the Evaluation part can be found in `./script/download.sh`.
You can download and generate the data sets and store them in a directory `integer_data`.

### Microbenchmark
The Microbenchmark of LeCo is located in `benchmark`, where you can setup by:
```
python3 fix_int_benchmark.py > fix_int.sh
bash fix_int.sh > fix_int_benchmark_intel.log
python3 auto_int_benchmark.py > auto_int.sh
bash auto_int.sh > auto_int_benchmark_intel.log
```
To parse the results, we provide the script in `scripts/paper_scripts/plot_cr_ra.ipynb`.
Except for system evaluation Figures, all other figures used in our paper is included in the above script.

### Layout
The file layout and explainations of our repository is as follows:
```
Learn-to-Compress
       \_________ benchmark 
              \____ fix_int (benchmarking FOR, Delta_fix, LeCo_fix)
              \____ auto_int (benchmarking Delta_var, LeCo_var)
       \_________ experiments (including all experiment of microbenchmark, as well as String Extensions)
       \_________ headers (implementation of different compression schemes)
       \_________ scripts (Setup/Data generation/Figure plot)
       \_________ thirdparty
              \____ succinct (Elias-Fano)
              \____ FSST (String Baseline)
       
```
