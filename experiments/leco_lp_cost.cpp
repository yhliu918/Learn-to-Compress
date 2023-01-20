#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"



int random(int m)
{
  return rand() % m;
}

void parse_segment_info(std::string filename, std::vector<int> &start_index, std::vector<double> &theta0, std::vector<double> &theta1)
{
  std::ifstream srcFile(filename, std::ios::in);
  while (srcFile.good())
  {
    int tmp;
    int tmp0;
    double tmp1;
    double tmp2;
    srcFile >> tmp >> tmp0 >> tmp >> tmp1 >> tmp2;
    start_index.push_back(tmp0);
    theta0.push_back(tmp1);
    theta1.push_back(tmp2);
    if (!srcFile.good())
    {
      break;
    }
  }
  srcFile.close();
}


int main(int argc, const char* argv[]) {
  using namespace Codecset;
  std::string method = std::string(argv[1]);
  std::string source_file = std::string(argv[2]);
  int delta = atoi(argv[3]);
  int model_size = atoi(argv[4]);

  // We pick a CODEC
  IntegerCODEC& codec = *CODECFactory::getFromName(method);

  std::vector<uint32_t> data;
  std::ifstream srcFile("../data/"+source_file+".txt", std::ios::in);
  if (!srcFile) {
    std::cout << "error opening source file." << std::endl;
    return 0;
  }
  int counter = 0;
  while (srcFile.good()) {
    counter++;
    uint32_t next;
    srcFile >> next;
    if (!srcFile.good()) { break; }
    data.push_back(next);

  }
  srcFile.close();
  int N = data.size();
  if (data.size() == 0) {
    std::cout << "Empty vector" << std::endl;
    return 0;
  }

  std::vector<int> start_index;
  std::vector<double> theta0;
  std::vector<double> theta1;
  parse_segment_info("../scripts/leco_lp_cost/"+source_file+".log", start_index, theta0, theta1);
  start_index.push_back(N);
  int blocks = start_index.size() - 1;
  int block_size = data.size() / blocks;


  codec.init(blocks, block_size, delta);
  int totalsize = 0;
  for (int i = 0;i < blocks;i++) {
    int block_length = start_index[i+1] - start_index[i];
    uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t));
    uint8_t *res = descriptor;
    res = codec.encodeArray8(data.data() + start_index[i], block_length, res, i);
    descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
    totalsize += (res - descriptor);
  }

  double origin_size = (sizeof(uint32_t) * N * 1.0);
  double compressrate = (totalsize)*100.0 / origin_size;
  std::cout << "compressed size " << totalsize << " uncompressed size " << 4 * N << std::endl;
  std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;

  codec.destroy();



}
