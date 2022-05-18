// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/piecewise_fix_string.h"
#include <gperftools/profiler.h>
//#include "headers/string/piecewise_fix_string_modify.h"
using namespace Codecset;


std::vector<uint8_t *> block_start_vec;
std::string firstkey_each_block;
int block_size;
int max_string;
int string_len;

bool pre_Bsearch(int left, int right, int *index, std::string &key, bool* skip_linear_scan)
{ // check which block the key belongs to
  while (left != right)
  {
    // The `mid` is computed by rounding up so it lands in (`left`, `right`].
    int64_t mid = left + (right - left + 1) / 2;
    std::string tmp_key = firstkey_each_block.substr(string_len*mid, string_len);
    if (tmp_key < key)
    {
      // Key at "mid" is smaller than "target". Therefore all
      // blocks before "mid" are uninteresting.
      left = mid;
    }
    else if (tmp_key > key)
    {
      // Key at "mid" is >= "target". Therefore all blocks at or
      // after "mid" are uninteresting.
      right = mid - 1;
    }
    else
    {
      *skip_linear_scan = true;
      left = right = mid;
    }
  }

  if (left == -1)
  {
    // All keys in the block were strictly greater than `target`. So the very
    // first key in the block is the final seek result.
    *skip_linear_scan = true;
    *index = 0;
  }
  else
  {
    *index = static_cast<uint32_t>(left);
  }
  return true;
}

bool Bsearch(int low, int high, int index,long_int &key)
{
  int mid;
  if (low > high)
  {
    return false;
  }
  mid = (low + high) / 2;
  mpz_t data_mid;
  mpz_init(data_mid);
  randomdecodeArray8_longint(block_start_vec[index], mid % block_size, &data_mid);
  long_int tmp_key(data_mid);
  if (tmp_key == key)
  {
    return true;
  }
  else if (tmp_key > key)
  {
    return Bsearch(low, mid - 1, index, key);
  }
  else if (tmp_key < key)
  {
    return Bsearch(mid + 1, high, index, key);
  }
}

bool ourBsearch(int low, int high, long_int &key)
{
  int mid;
  if (low > high)
  {
    return false;
  }
  mid = (low + high) / 2;
  mpz_t data_mid;
  mpz_init(data_mid);
  randomdecodeArray8_longint(block_start_vec[(int)mid / block_size], mid % block_size, &data_mid);
  long_int tmp_key(data_mid);
  if (tmp_key == key)
  {
    return true;
  }
  else if (tmp_key > key)
  {
    return ourBsearch(low, mid - 1, key);
  }
  else if (tmp_key < key)
  {
    return ourBsearch(mid + 1, high, key);
  }
  return false;
}

int main()
{

  // We pick a CODEC

  // could use others, e.g., "simdfastpfor256", "BP32"
  ProfilerStart("test_capture.prof");
  std::vector<std::string> string_vec;
  std::ifstream srcFile("/home/lyh/rocksdb/dump_data/padding_a_wholestring_key.txt", std::ios::in);
  // std::ifstream srcFile("/home/lyh/Learn-to-Compress/data/poisson_2000000/key.txt", std::ios::in);
  if (!srcFile)
  {
    std::cout << "error opening source file." << std::endl;
    return 0;
  }
  int cnt = 0;
  while (1)
  {
    std::string tmp_str;
    srcFile >> tmp_str;
    if (srcFile.eof())
    {
      break;
    }
    // std::cout << next << std::endl;
    string_vec.push_back(tmp_str);
  }

  srcFile.close();

  int string_length = string_vec[0].size();
  int N = string_vec.size();
  int blocks = 100;
  block_size = N / blocks;
  if (block_size * blocks < N)
  {
    blocks++;
  }
  //=========================
  // long_int tmp_val = convertToASCII(string_vec[0]);
  // std::cout<<string_vec[0]<<" "<<tmp_val<<" "<<convertToString(tmp_val)<<std::endl;
  //=========================

  std::vector<int> start_index;
  int totalsize = 0;
  max_string = 8 * string_vec[0].size();
  string_len = string_vec[0].size();
  std::vector<long_int> data_vec;
  for (int i = 0; i < blocks; i++)
  {
    int block_length = block_size;
    if (i == blocks - 1)
    {
      block_length = N - (blocks - 1) * block_size;
    }

    uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t) * 1000);
    uint8_t *res = descriptor;

    res = encodeArray8_string(string_vec, i * block_size, block_length, descriptor, max_string);
    // std::cout<<"block "<<i<<" "<<res-descriptor<<std::endl;
    descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
    block_start_vec.push_back(descriptor);
    totalsize += (res - descriptor);
  }

  for (int i = 0; i < blocks; i++)
  {
    firstkey_each_block.append(string_vec[i * block_size]);
  }

  double compressrate = (totalsize)*100.0 / (N * string_length * 1.0);
  std::cout << N << " " << string_length << " " << totalsize << " " << compressrate << std::endl;
  std::cout << "compressed size " << totalsize << " uncompressed size " << N * string_length << std::endl;
  std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;
  // long_int record = convertToASCII(string_vec[0]);
  // std::cout<< convertToString(record)<<std::endl;

  bool flag = true;

  double totaltime = 0.0;
  std::cout << "decompress all!" << std::endl;
  double start = getNow();
  for (int i = 0; i < blocks; i++)
  {
    int block_length = block_size;
    if (i == blocks - 1)
    {
      block_length = N - (blocks - 1) * block_size;
    }

    std::vector<std::string> result;

    decodeArray8_string(block_start_vec[i], block_length, max_string, result);

    // for (int j = 0; j < block_length; j++)
    // {

    //   if (string_vec[j + block_size * i].compare(result[j]))
    //   {
    //     std::cout << "num: " << j << " true is: " << string_vec[j + block_size * i] << " predict is: " << result[j] << std::endl;
    //     std::cout << "something wrong! decompress failed" << std::endl;
    //     flag = false;
    //     break;
    //   }
    // }
  }

  double end = getNow();
  totaltime += (end - start);

  std::cout << "all decoding time per int: " << std::setprecision(8)
            << totaltime / N * 1000000000 << " ns" << std::endl;

  std::cout << "random access decompress!" << std::endl;
  double randomaccesstime = 0.0;
  start = getNow();

  for (int i = 0; i < N; i++)
  {
    int index = i;
    std::string tmpvalue;
    randomdecodeArray8(block_start_vec[(int)index / block_size], index % block_size, tmpvalue);
    //std::cout<<i<<" "<<tmpvalue<<" block "<<(int)index / block_size<<std::endl;
     if (string_vec[index].compare(tmpvalue))
     {
       std::cout << "num: " << index << "true is: " << string_vec[index] << " predict is: " << tmpvalue << std::endl;
       flag = false;
       std::cout << "something wrong! decompress failed" << std::endl;
     }
     if (!flag)
     {
       break;
     }
  }
  end = getNow();
  randomaccesstime += (end - start);

  std::cout << "random decoding time per int: " << std::setprecision(8)
            << randomaccesstime / N * 1000000000 << " ns" << std::endl;

  int sample_size = 1000;
  start = getNow();
  for (int i = 0; i < sample_size; i++)
  {
    int index = random(N);
    std::string tmpkey = string_vec[index];
    bool skip_linear = false;
    int index_search = 0;
    pre_Bsearch(-1, blocks -1,&index_search,tmpkey,&skip_linear);
    //std::cout<< "index_search "<<index_search<<std::endl;
    if(!skip_linear){
        long_int record = convertToLongInt(tmpkey);
        int lower = index_search * block_size;
        int upper = (index_search + 1) * block_size;
        if(index_search == blocks -1){
          upper = N;
        }
        Bsearch(lower, upper, index_search,record);
    }
    //ourBsearch(0, N, tmp_val);
  }
  end = getNow();
  double ourbinarytime = end - start;
  std::cout << "binary time per time: " << std::setprecision(8)
            << ourbinarytime / sample_size * 1000000000 << " ns" << std::endl;

  ProfilerStop();
}
