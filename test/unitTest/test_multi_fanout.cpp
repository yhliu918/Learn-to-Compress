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
#include "../headers/piecewise_multi_fanout.h"


int main() {
  using namespace Codecset;

  // We pick a CODEC
  
  IntegerCODEC &codec = *CODECFactory::getFromName("piecewise_multi_fanout");
  (dynamic_cast<piecewise_multi_fanout*>(&codec))->seg_num = 4; //can be customized

  std::vector<uint32_t> data;

  std::ifstream srcFile("../data/standard/linear_200M_uint32.txt", std::ios::in);
  if(!srcFile) { 
      std::cout << "error opening source file." << std::endl;
      return 0;
  }
  while(1){     
      uint32_t next ;
      srcFile  >> next;
      if(srcFile.eof()){break;}
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
  
  int blocks = 1000;
  int block_size = data.size()/blocks; 
  int delta  = 1; 
    
  //std::vector<uint8_t*> block_start_vec;
  //int totalsize = 0;
  codec.init(blocks, block_size, delta);


  for(int i=0;i<blocks;i++){
    //uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t));
    //uint8_t * res = descriptor;
    uint8_t * tmp_pointer = NULL;
    codec.encodeArray8(data.data()+(i*block_size),block_size ,tmp_pointer,i);
    //descriptor = (uint8_t*)realloc(descriptor, (res-descriptor));
    //block_start_vec.push_back(descriptor);
    //totalsize += (res-descriptor); 
  }
  
  double compressrate = ((dynamic_cast<piecewise_multi_fanout*>(&codec))->total_byte) *100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;

  bool flag =true;
  std::vector<uint32_t> recover(data.size());
  double totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
  double start = getNow();
  uint8_t *tmp_pointer = NULL;
  codec.decodeArray8(tmp_pointer, block_size, recover.data(), 0);
  for(int i = 0; i < N; ++i) {
    if(recover[i] != data[i]) {
      std::cout << "no. " <<i<< "true is: "<<data[i]<<" decode result is " << recover[i] <<std::endl;
      std::cout << "something wrong! decode failed"<<std::endl;
      flag = false;
    }
    if(!flag) break;
  }
  double end = getNow();
  totaltime += (end - start);

  std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
  std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;
  
  std::vector<uint32_t> buffer(data.size());
  double randomaccesstime =0.0;
  std::cout<<"random access decompress!"<<std::endl; 
  start = getNow();
  uint32_t mark=0;
    
  for(int i=0;i<N;i++){ 
      uint8_t *tmp_pointer =NULL;  
      uint32_t tmpvalue = codec.randomdecodeArray8(tmp_pointer, i, buffer.data(), 0);
       mark+=tmpvalue;

       if(data[i]!=tmpvalue){
        
        std::cout<<"num: "<<i<< "true is: "<<data[i]<<" predict is: "<<tmpvalue<<std::endl;
        flag = false;
        std::cout<<"something wrong! decompress failed"<<std::endl;
        
      }
    if(!flag){
        break;
    }
    }
  end = getNow();
  randomaccesstime+=(end-start);
  
  std::cout << "random decoding time per int: " << std::setprecision(8)
     << randomaccesstime / N * 1000000000 << "ns" << std::endl;
  std::cout << "random decoding speed: " << std::setprecision(10)
     << N /(randomaccesstime*1000) <<  std::endl;

}
