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
#include "../headers/split.h"






int main() {
  using namespace Codecset;

  // We pick a CODEC
  IntegerCODEC &codec = *CODECFactory::getFromName("piecewise_fanout");

  std::vector<uint32_t> data;
  std::ifstream srcFile("/root/Learn-to-Compress/data/wf/wiki.txt",std::ios::in); 
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
 
    
  int block_size = 6300;
  int blocks = data.size()/block_size;
  if(block_size*blocks<N){
    blocks++;
  }
  codec.init(blocks,block_size,0);
  int totalsize = 0;
  uint8_t * res = NULL;

  for(int i=0;i<blocks;i++){
    res = codec.encodeArray8(data.data()+(i*block_size),block_size ,res,i);
  }
  totalsize = codec.get_block_nums();
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;
  
  
  bool flag =true;
  uint32_t * recover = new uint32_t[data.size()];
  double totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
   double start = getNow();
   
    codec.decodeArray8(res, block_size, recover, N);
      
      for(int j=0;j<N;j++){
        if(data[j]!=recover[j]){
          std::cout<<"num: "<<j<< " true is: "<<data[j]<<" predict is: "<<recover[j]<<std::endl;
          std::cout<<"something wrong! decompress failed"<<std::endl;
          flag = false;
          break;
         }
         
       }
       
       


      double end = getNow();
      totaltime += (end - start);

std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;

    

  std::cout<<"random access decompress!"<<std::endl; 
  
  double randomaccesstime =0.0;
   start = getNow();
   uint32_t mark=0;
   uint32_t* placeholder=NULL;
    
  for(int i=0;i<N;i++){    
      
      uint32_t tmpvalue = codec.randomdecodeArray8(res, i,placeholder , i/block_size);
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
     << randomaccesstime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "random decoding speed: " << std::setprecision(10)
     << data.size()/(randomaccesstime*1000) <<  std::endl;
    
  codec.destroy();  



  
}
