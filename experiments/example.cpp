// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/piecewise_fix_string_outlier_detect.h"
//#include "headers/string/piecewise_fix_string_modify.h"
int random(int m)
{
  return rand() % m;
}

int main()
{
  using namespace Codecset;

  // We pick a CODEC

  // could use others, e.g., "simdfastpfor256", "BP32"

  std::vector<std::string> string_vec = {"abcdefgzaaa", "abcdefgzacd", "abcdefgzaef", "abcdefgzagh","abcdefgzaij","abcdefgzakl", "abcdefgzzzz"};
  int string_length = string_vec[0].size();
  int N = string_vec.size();
  int blocks = 1;
  int block_size = N / blocks;

  std::vector<uint8_t *> block_start_vec;
  std::vector<int> start_index;
  int totalsize = 0;
  for (int i = 0; i < blocks; i++)
  {
    int block_length = block_size;
    if (i == blocks - 1)
    {
      block_length = N - (blocks - 1) * block_size;
    }

    uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t) * 100);
    uint8_t *res = descriptor;
    res = encodeArray8_string(string_vec, block_length, descriptor, 11*8);

    descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
    block_start_vec.push_back(descriptor);
    totalsize += (res - descriptor);
  }

  double compressrate = (totalsize)*100.0 / (N * string_length * 1.0);
  std::cout << "total compression rate:" << std::setprecision(4) << compressrate << std::endl;
  // long_int record = convertToASCII(string_vec[0]);
  // std::cout<< convertToString(record)<<std::endl;

  bool flag = true;
  long_int *recover = new long_int[N];

  double totaltime = 0.0;
  std::cout << "decompress all!" << std::endl;
  double start = getNow();
  for(int i=0;i<blocks;i++){
    std::vector<std::string> result;
    decodeArray8_string(block_start_vec[i], block_size, recover, 11*8,result);
    for (int j = 0; j < N; j++)
    {

      if (string_vec[j].compare(result[j]))
      {
        std::cout << "num: " << j << " true is: " << string_vec[j] << " predict is: " << result[j] << std::endl;
        std::cout << "something wrong! decompress failed" << std::endl;
        flag = false;
        break;
      }
    }


  }

  double end = getNow();
  totaltime += (end - start);

  std::cout << "all decoding time per int: " << std::setprecision(8)
            << totaltime /N * 1000000000 << "ns" << std::endl;

}
