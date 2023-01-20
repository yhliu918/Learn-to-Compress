#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"



int random(int m)
{
  return rand() % m;
}


int main(int argc, const char* argv[]) {
  using namespace Codecset;
  std::string method = std::string(argv[1]);
  std::string source_file = std::string(argv[2]);
  int blocks = atoi(argv[3]);
  int model_size = atoi(argv[4]);

  // We pick a CODEC
  IntegerCODEC& codec = *CODECFactory::getFromName(method);

  std::vector<uint32_t> data;
  std::ifstream srcFile("../data/"+source_file, std::ios::in);
  if (!srcFile) {
    std::cout << "error opening source file." << std::endl;
    return 0;
  }
  int counter = 0;
  int cut = 0;
  while (srcFile.good()) {
    
    counter++;
    uint32_t next;
    srcFile >> next;
    if (!srcFile.good()) { break; }
    if(counter>cut){
      data.push_back(next);
    } 

  }
  srcFile.close();
  int N = data.size();
  if (data.size() == 0) {
    std::cout << "Empty vector" << std::endl;
    return 0;
  }

  
  int blocks_real = N;
  int block_size = data.size() / blocks;
  // int delta = (1<<30);
  codec.init(blocks, block_size, 0);
  int totalsize = 0;
  uint8_t* res = NULL;
  int total_seg = 0;
  int block_length = block_size;
  for (int i = 0;i < blocks;i++) {
    // std::cout<<i<<std::endl;
    if(i==blocks-1){
      block_length = N - i*block_size;
    }
    res = codec.encodeArray8(data.data() + (i * block_size), block_length, res, i);
    
  }
  totalsize = codec.get_block_nums();

  
  double origin_size = (sizeof(uint32_t) * N * 1.0);
  double compressrate = (totalsize)*100.0 / origin_size;
  std::cout << "compressed size " << totalsize << " uncompressed size " << 4 * N << std::endl;
  std::cout << "total compression rate: " << std::setprecision(4) << compressrate << std::endl;
  bool flag = true;

  std::vector<uint32_t> recover(data.size());
  double totaltime = 0.0;
  
  // std::cout<<"decompress all!"<<std::endl;
  //  double start = getNow();
  // for(int i=0;i<blocks;i++){

  //     codec.decodeArray8(res, block_size, recover.data()+i*block_size, i);

  //     for(int j=0;j<block_size;j++){
  //       if(data[j+i*block_size]!=recover[j+i*block_size]){
  //         std::cout<<"block: "<<i<<" num: "<<j<< " true is: "<<data[j+i*block_size]<<" predict is: "<<recover[j+i*block_size]<<std::endl;
  //         std::cout<<"something wrong! decompress failed"<<std::endl;
  //         flag = false;
  //         break;
  //        }

  //      }
  //      if(!flag){
  //         break;
  //      }


  // }
  //     double end = getNow();
  //     totaltime += (end - start);
  // double da_ns = totaltime / data.size() * 1000000000;

// std::cout << "all decoding time per int: " << std::setprecision(8)
//      << totaltime / data.size() * 1000000000 << "ns" << std::endl;
// std::cout << "all decoding speed: " << std::setprecision(10)
//      << data.size()/(totaltime*1000) <<  std::endl;

  // std::cout << "random access decompress!" << std::endl;
  std::vector<uint32_t> ra_pos;
  for(int i=0;i<N;i++){
    ra_pos.push_back(random(N));
  }
  double randomaccesstime = 0.0;
  double start = getNow();
  uint32_t mark = 0;
  uint32_t* placeholder = NULL;

  for (auto index: ra_pos) {

    uint32_t tmpvalue = codec.randomdecodeArray8(res, index, placeholder, index / block_size);
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
  double ra_ns = randomaccesstime / N * 1000000000;

  std::cout<<method<<" "<<source_file<<" "<<blocks<<" "<<compressrate<<" "<<ra_ns<<std::endl;
  // std::cout << "random decoding time per int: " << std::setprecision(8)
  //   << randomaccesstime / data.size() * 1000000000 << " ns" << std::endl;
  // std::cout << "random decoding speed: " << std::setprecision(10)
  //   << data.size() / (randomaccesstime * 1000) << std::endl;

  codec.destroy();
  // free(&codec);



}
