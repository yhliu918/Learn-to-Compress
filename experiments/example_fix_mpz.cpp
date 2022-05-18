// searching for hyper-parameter like block number
#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/string_utils.h"
#include "../headers/string/piecewise_fix_string_mpz_t.h"

using namespace Codecset;

int main()
{

  // We pick a CODEC

  // could use others, e.g., "simdfastpfor256", "BP32"
  std::vector<std::string> string_vec;
  // std::ifstream srcFile("/home/lyh/string_data/email_list/padding_0_prefix.txt", std::ios::in);
  std::ifstream srcFile("/home/lyh/Learn-to-Compress/data/poisson_20000/key.txt", std::ios::in);
  if (!srcFile)
  {
    std::cout << "error opening source file." << std::endl;
    return 0;
  }
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

  std::vector<uint8_t *> block_start_vec;

  int string_length = string_vec[0].size();
  int N = string_vec.size();
  int blocks = 100;
  int block_size = N / blocks;
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
  int max_string = 8 * string_vec[0].size();
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

    descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
    block_start_vec.push_back(descriptor);
    totalsize += (res - descriptor);
  }

  double compressrate = (totalsize)*100.0 / (N * string_length * 1.0);
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
}
