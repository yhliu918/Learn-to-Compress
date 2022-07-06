import os
import glob
from pathlib import Path
import datetime

method_list = ["FOR", "delta_my", "piecewise_fix_op"]
source_file_and_blocks = {
    "books_200M_uint32.txt": [1000000, 1000000, 1000000],
    "normal_200M_uint32.txt": [1000000, 1000000, 1000000],
    "linear_200M_uint32.txt": [1000000, 1000000, 1000000],
    "fb/fb-289000.txt": [1000, 1000, 1000],
    "wf/wiki.txt": [10000, 10000, 10000],
    "wf/newman.txt": [1000, 1000, 1000]
}
bitmap_names_200M = ["bitmap_random_1e-05_200000000.txt", "bitmap_random_0.0001_200000000.txt", "bitmap_random_0.0005_200000000.txt",
                     "bitmap_random_0.001_200000000.txt", "bitmap_random_0.01_200000000.txt", "bitmap_random_0.1_200000000.txt"]
bitmap_names_289K = ["bitmap_random_1e-05_289000.txt", "bitmap_random_0.0001_289000.txt", "bitmap_random_0.0005_289000.txt",
                     "bitmap_random_0.001_289000.txt", "bitmap_random_0.01_289000.txt", "bitmap_random_0.1_289000.txt"]
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
output = Path(timestamp)
os.mkdir(output)


def run_exp():
    curdir = os.getcwd()
    print(curdir)
    if curdir.split("/")[-1] != "Learn-to-Compress":
        exit(1, "Please run this script in Learn-to-Compress folder")
    os.system("cmake --build ./build --config Release --target test_bitmap --")

    for source_file in source_file_and_blocks:
        blocks_num_list = source_file_and_blocks[source_file]
        for i, method in enumerate(method_list):
            b_list = bitmap_names_289K if source_file == "fb/fb-289000.txt" else bitmap_names_200M
            for bitmap in b_list:
                os.system(
                    f"./build/test_bitmap {method} {source_file} {bitmap} {str(blocks_num_list[i])} >> {str(output)}/exp_{method}_{source_file.split('/')[-1]}.txt")
            print("Done with " + method + " for " + source_file)
        print("Done with all methods for " + source_file)


def post_process(out_dir):
    result_file_list = glob.glob(f"{out_dir}/exp_*.txt")
    for result_file in result_file_list:
        print(result_file)
        with open(result_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("access time per int"):
                    print(line.split(":")[-1].strip().split(" ")[0])


# main
if __name__ == "__main__":
    run_exp()
    # post_process('20220706_171206')
    post_process(str(output))
