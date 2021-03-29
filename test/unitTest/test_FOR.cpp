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






int main() {
  using namespace Codecset;

  // We pick a CODEC
  IntegerCODEC &codec = *CODECFactory::getFromName("FOR");
  
  // could use others, e.g., "simdfastpfor256", "BP32"


  std::vector<uint32_t> data;
  std::ifstream srcFile("../data/fb/fb-289000.txt",std::ios::in); 
  if(!srcFile) { 
      std::cout << "error opening source file." << std::endl;
      return 0;
  }
  while(1){
      
      uint32_t next ;
      srcFile >> next;
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

  int blocks=1000;
  int delta =0;
  int block_size = data.size()/blocks;
  bool flag = true;
  codec.init(blocks,block_size,delta);
  std::vector<uint8_t*> block_start_vec;
  int totalsize=0;
  for(int i=0;i<blocks;i++){
    uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t)+1024);
    uint8_t * res = descriptor;
    res = codec.encodeArray8(data.data()+i*block_size,block_size, descriptor,N);
    descriptor = (uint8_t*)realloc(descriptor, (res-descriptor)*sizeof(uint8_t));
    block_start_vec.push_back(descriptor);
    totalsize += (res-descriptor);
  }

  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;
    
  std::vector<uint32_t> recover(data.size());
  double totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
  double start = getNow();
  for(int i=0;i<blocks;i++){
      codec.decodeArray8(block_start_vec[i], block_size, recover.data()+i*block_size, i);
      /*
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
       */
  }
double end = getNow();
totaltime += (end - start);

std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;

  std::cout<<"random access decompress!"<<std::endl; 
  std::vector<uint32_t> buffer(data.size());
  double randomaccesstime =0.0;
  start = getNow();
  for(int i=0;i<N;i++){    

      uint32_t tmpvalue = codec.randomdecodeArray8(block_start_vec[i/block_size], i%block_size, buffer.data(), i/block_size);
      

      //std::cout<<"processing...  "<<j<<" / "<<N<<std::endl;
       /*
      if(data[i]!=tmpvalue){
        
        std::cout<<"num: "<<i<< "true is: "<<data[i]<<" predict is: "<<tmpvalue<<std::endl;
        flag = false;
        std::cout<<"something wrong! decompress failed"<<std::endl;
        
      }
    if(!flag){
        break;
    }
     */
  }
end = getNow();
randomaccesstime+=(end-start);


std::cout << "random decoding time per int: " << std::setprecision(8)
     << randomaccesstime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "random decoding speed: " << std::setprecision(10)
     << data.size()/(randomaccesstime*1000) <<  std::endl;





  
}
