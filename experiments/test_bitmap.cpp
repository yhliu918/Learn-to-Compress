
#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
#include "../headers/bitvector.hpp"

int main(int argc, const char *argv[])
{
  using namespace Codecset;
  std::string method = std::string(argv[1]);
  std::string source_file = std::string(argv[2]);
  std::string bitmap_name = std::string(argv[3]);
  int blocks = atoi(argv[4]);
  // We pick a CODEC
  // IntegerCODEC &codec = *CODECFactory::getFromName("delta_my");
  // IntegerCODEC &codec = *CODECFactory::getFromName("piecewise_fix_op");
  IntegerCODEC &codec = *CODECFactory::getFromName(method);

  std::vector<uint32_t> data;
  std::ifstream srcFile("/home/lyh/Learn-to-Compress/integer_data/" + source_file, std::ios::in);
  if (!srcFile)
  {
    std::cout << "error opening source file." << std::endl;
    return 0;
  }
  while (1)
  {

    uint32_t next;
    srcFile >> next;
    if (srcFile.eof())
    {
      break;
    }
    data.push_back(next);
  }
  srcFile.close();
  int N = data.size();
  std::vector<uint32_t> bitmap;
  std::ifstream bitFile("/home/zxy/Learn-to-Compress-0510/data/bitmap_random/" + bitmap_name, std::ios::in);
  std::cout << "../data/bitmap_random/" + bitmap_name << std::endl;
  int k = 0;
  std::vector<uint32_t> bit_pos;
  for (int i = 0; i < N; i++)
  {

    uint32_t next;
    bitFile >> next;
    if (next)
    {
      bit_pos.emplace_back(i);
    }
    k++;
    if (bitFile.eof())
    {
      break;
    }
    bitmap.push_back(next);
  }
  bitFile.close();
  Bitvector bitvector(bitmap);

  int block_size = data.size() / blocks;
  blocks = data.size() / block_size;
  if (blocks * block_size < N)
  {
    blocks++;
  } // handle with the last block, maybe < block_size
  std::cout << "Total blocks " << blocks << " block size " << block_size << std::endl;
  int delta = 32;
  codec.init(blocks, block_size, delta);
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
    uint8_t *descriptor = (uint8_t *)malloc(block_length * sizeof(uint64_t));
    uint8_t *res = descriptor;
    res = codec.encodeArray8(data.data() + (i * block_size), block_length, descriptor, N);
    descriptor = (uint8_t *)realloc(descriptor, (res - descriptor));
    block_start_vec.push_back(descriptor);
    totalsize += (res - descriptor);
  }

  double compressrate = (totalsize)*100.0 / (4 * N * 1.0);
  std::cout << "total compression rate:" << std::setprecision(4) << compressrate << std::endl;
  //************************************************************
  /*
    uint32_t* buffer = NULL;
    double start = getNow();
    uint32_t tmpvalue =0;
    int counter = 0;
    for(int i=0;i<N;i++){
      if(bitmap[i]){
        counter++;
        tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)i/block_size], i%block_size, buffer, i/block_size);
      }
    }
    double end = getNow();
    std::cout<<"counter "<<counter<<std::endl;
    std::cout << "access time per int: " << std::setprecision(8)
       << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;
  //************************************************************
    for(int i=0;i<(int)block_start_vec.size();i++){
         free(block_start_vec[i]);
     }
    }
    */

  // PIECEWISE
  uint32_t *buffer = NULL;
  double start = getNow();
  uint32_t tmpvalue = 0;
  int counter = 0;
  for (auto &pos : bit_pos)
  {
    counter++;
    tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)pos / block_size], pos % block_size, buffer, N);
  }
  double end = getNow();
  std::cout << bit_pos.size() << std::endl;
  std::cout << "access time per int: " << std::setprecision(8)
            << (end - start) / N * 1000000000 << " ns" << std::endl;

  /*
  // PIECEWISE DETECT
    uint32_t* buffer = NULL;
    std::vector<uint32_t> recover(block_size);
    double start = getNow();
    uint32_t tmpvalue =0;
    double  total_count = 0;
    for(int i=0;i<blocks;i++){
        //std::cout<<i<<std::endl;
        int total_one =0;
        for(int j=0;j<block_size;j++){
            if(bitmap[i*block_size+j]){total_one++;}
        }
        if(total_one>block_size*0.198){
            total_count++;
            codec.decodeArray8(block_start_vec[i], block_size, recover.data(), i);
            for(int j=0;j<block_size;j++){
              if(bitmap[i*block_size+j]){
                tmpvalue = recover[j];
              }
            }
        }
        else{
          for(int j=0;j<block_size;j++){
              if(bitmap[i*block_size+j]){
                tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)i], j, buffer, i);
              }
            }
        }
    }
    double end = getNow();
    std::cout << "access time per int: " << std::setprecision(8)
       << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;
    std::cout<< "about "<< total_count<<" blocks use decode all" <<std::endl;
  */
  /*
    // DELTA DETECT
    std::vector<uint32_t> recover(block_size);
    double start = getNow();
    uint32_t tmpvalue =0;
    for(int i=0;i<blocks;i++){
        //std::cout<<i<<std::endl;
        bool flag = false;
        for(int j=0;j<block_size;j++){
          if(bitmap[i*block_size+j]){
            if(flag){
              tmpvalue = recover[j];
            }
            else{
              codec.decodeArray8(block_start_vec[i], block_size, recover.data(), i);
              tmpvalue = recover[j];
              flag = true;
            }
          }
        }
    }
    double end = getNow();
    std::cout << "access time per int: " << std::setprecision(8)
       << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;
  */

  /*
    // DELTA DIRECT
    std::vector<uint32_t> recover(block_size);
    double start = getNow();
    uint32_t tmpvalue =0;
    for(int i=0;i<blocks;i++){
        codec.decodeArray8(block_start_vec[i], block_size, recover.data(), i);
        for(int j=0;j<block_size;j++){
          if(bitmap[i*block_size+j]){
              tmpvalue = recover[j];
          }
        }
    }
    double end = getNow();
    std::cout << "access time per int: " << std::setprecision(8)
       << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;
  */
  for (int i = 0; i < (int)block_start_vec.size(); i++)
  {
    free(block_start_vec[i]);
  }
}
