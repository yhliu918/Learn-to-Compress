// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/string/lr_string.h"
#include "../headers/string/string_utils.h"
//#include "../headers/string/piecewise_fix_string.h"
#include "../headers/string/piecewise_auto_string.h"


int main()
{
  using namespace Codecset;

  // We pick a CODEC

  // could use others, e.g., "simdfastpfor256", "BP32"

  std::vector<std::string> string_vec;
  std::ifstream srcFile("/home/lyh/Learn-to-Compress/data/poisson_20000/key.txt", std::ios::in);
  // std::ifstream srcFile("/home/zxy/Learn-to-Compress/data/result_key.txt", std::ios::in);
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

  int string_length = string_vec[0].size();
  int N = string_vec.size();
  long_int delta = ((long_int)1 << 35);
  Piecewise_auto codec;
  codec.init(delta);

  int max_string = 8 * string_vec[0].size();
  codec.encodeArray8(string_vec, 0, N, NULL, max_string);
  int totalsize = codec.get_total_byte();

  double compressrate = (totalsize)*100.0 / (N * string_length * 1.0);
  std::cout<<"compressed size "<<totalsize<<" uncompressed size "<<N * string_length<<std::endl;
  std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;
  std::cout << "total segment number: " << codec.get_total_seg() << std::endl;
  // long_int record = convertToASCII(string_vec[0]);
  // std::cout<< convertToString(record)<<std::endl;

  bool flag = true;

  double totaltime = 0.0;
  std::cout << "decompress all!" << std::endl;
  double start = getNow();

  std::vector<std::string> result;
  long_int *recover = new long_int[N];
  codec.decodeArray8(NULL, N, recover, max_string, result);

  // for (int j = 0; j < N; j++)
  // {

  //   if (string_vec[j].compare(result[j]))
  //   {
  //     std::cout << "num: " << j << " true is: " << string_vec[j] << " predict is: " << result[j] << std::endl;
  //     std::cout << "something wrong! decompress failed" << std::endl;
  //     flag = false;
  //     break;
  //   }
  // }

  double end = getNow();
  totaltime += (end - start);

  std::cout << "all decoding time per int: " << std::setprecision(8)
            << totaltime / N * 1000000000 << "ns" << std::endl;
}
