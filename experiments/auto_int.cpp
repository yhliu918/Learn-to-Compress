//
// A simple example to get you started with the library.
// You can compile and run this example like so:
//
//   make example
//   ./example
//
//  Warning: If your compiler does not fully support C++11, some of
//  this example may require changes.
//

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"



int random(int m)
{
  return rand() % m;
}


int main() {
  using namespace Codecset;

  // We pick a CODEC
  IntegerCODEC& codec = *CODECFactory::getFromName("piecewise_cost_dp");

  std::vector<uint32_t> data;
  std::ifstream srcFile("/home/lyh/Learn-to-Compress/integer_data/fb/fb-289000.txt", std::ios::in);
  if (!srcFile) {
    std::cout << "error opening source file." << std::endl;
    return 0;
  }
  int counter = 0;
  int cut = 3000;
  while (srcFile.good()) {
    if(counter==cut){
      break;
    } 
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
  std::cout << "vector size = " << data.size() << std::endl;
  std::cout << "vector size = " << data.size() * sizeof(uint32_t) / 1024.0 << "KB"
    << std::endl;


  int blocks = 1;
  int blocks_real = N;
  int block_size = data.size() / blocks;
  int delta = (1<<11);
  // int delta = (1<<30);
  codec.init(blocks, block_size, delta);
  int totalsize = 0;
  uint8_t* res = NULL;

  for (int i = 0;i < blocks;i++) {
    res = codec.encodeArray8(data.data() + (i * block_size), block_size, res, blocks_real);
  }
  totalsize = codec.get_block_nums();
  double compressrate = (totalsize) * 100.0 / (4 * N * 1.0);
  std::cout << "compressed size " << totalsize << " uncompressed size " << 4 * N << std::endl;
  std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;
  bool flag = true;
  uint32_t* recover = new uint32_t[data.size()];

  //std::vector<uint32_t> recover(data.size());
  double totaltime = 0.0;
  /*
  std::cout<<"decompress all!"<<std::endl;
   double start = getNow();
  for(int i=0;i<blocks;i++){

      codec.decodeArray8(res, block_size, recover+i*block_size, i);

      for(int j=0;j<block_size;j++){
        if(data[j+i*block_size]!=recover[j+i*block_size]){
          std::cout<<"block: "<<i<<" num: "<<j<< " true is: "<<data[j+i*block_size]<<" predict is: "<<recover[j+i*block_size]<<std::endl;
          std::cout<<"something wrong! decompress failed"<<std::endl;
          flag = false;
          break;
         }

       }
       if(!flag){
          break;
       }


  }
      double end = getNow();
      totaltime += (end - start);

std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;
*/
  std::cout << "random access decompress!" << std::endl;

  double randomaccesstime = 0.0;
  double start = getNow();
  uint32_t mark = 0;
  uint32_t* placeholder = NULL;

  for (int i = 0;i < N;i++) {
    int index = random(N);
    // int index = i; 


    uint32_t tmpvalue = codec.randomdecodeArray8(res, index % block_size, placeholder, index / block_size);
    mark += tmpvalue;

    if(data[index]!=tmpvalue){

    std::cout<<"num: "<<index<< "true is: "<<data[index]<<" predict is: "<<tmpvalue<<std::endl;
    flag = false;
    std::cout<<"something wrong! decompress failed"<<std::endl;

    }
    if(!flag){
      break;
    }

  }
  double end = getNow();
  randomaccesstime += (end - start);

  std::cout << "random decoding time per int: " << std::setprecision(8)
    << randomaccesstime / data.size() * 1000000000 << " ns" << std::endl;
  std::cout << "random decoding speed: " << std::setprecision(10)
    << data.size() / (randomaccesstime * 1000) << std::endl;

  codec.destroy();




}
