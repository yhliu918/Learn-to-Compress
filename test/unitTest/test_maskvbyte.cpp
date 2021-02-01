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

#include "./headers/common.h"
#include "./headers/codecfactory.h"
#include "./headers/caltime.h"






int main() {
  using namespace Codecset;

  // We pick a CODEC
  IntegerCODEC &codec = *CODECFactory::getFromName("MaskVByte");
  // could use others, e.g., "simdfastpfor256", "BP32"
  ////////////
  //
  // create a container with some integers in it
  //
  // for this example, we will not assume that the
  // integers are in sorted order
  //
  // (Note: You don't need to use a vector.)
  //

  int N = 200000000;
  int blocks =1;
  std::vector<uint32_t> data(N+1024);
  FILE *fpRead=fopen("../data/books_200M_uint32.txt","r");

  for(int i=0;i<N;i++)
   {
      fscanf(fpRead,"%u",&data[i]);
   }
  fclose(fpRead);
  std::cout << "[standard benchmark]" << std::endl;


  if (data.size() == 0) {
    std::cout << "Empty vector" << std::endl;
    return 0;
  }
  std::cout << "vector size = " << data.size() << std::endl;
  std::cout << "vector size = " << data.size() * sizeof(uint32_t) / 1024.0 << "KB"
       << std::endl;
  

  int block_size = data.size()/blocks;
  double compressrate=0;
  double totaltime =0.0;

  for(int i=0;i<blocks;i++){
    std::vector<uint8_t> compdata(8*data.size()/blocks ); 
    std::vector<uint32_t> buffer(data.size()/blocks +1024);
    size_t compsize = compdata.size();

    uint8_t *out = codec.encodeArray8(data.data(),N, compdata.data(),compsize);
 
    compressrate += (out - compdata.data())*100.0  / (4*data.size()*1.0/blocks);


    size_t recoveredsize = buffer.size();
    //std::cout<<"decompress all!"<<std::endl;
    double start = getNow();
    codec.decodeArray8(compdata.data(), N, buffer.data(), recoveredsize);
    double end = getNow();
    totaltime += (end - start);
    for(int j=0;j<block_size;j++){
      if(data[j+i*block_size]!=buffer[j]){
        std::cout<<"block: "<<i<<"num: "<<j<< "true is: "<<data[j]<<" predict is: "<<buffer[j]<<std::endl;
        std::cout<<"something wrong! decompress failed"<<std::endl;
        break;
      }

    }
    
    //std::cout<<"random access decompress!"<<std::endl;
    /*
    for(int j=0;j<block_size;j++){
      start = getNow();
      uint32_t tmp = i/block_size;
      uint32_t tmpvalue = codec.randomdecodeArray(compdata.data(), j, buffer.data(), recoveredsize);
      end = getNow();
      randomaccesstime+=(end-start);
      if(data[j+i*block_size]!=tmpvalue){
        
        std::cout<<"block: "<<i<<"num: "<<j<< "true is: "<<data[j+i*block_size]<<" predict is: "<<tmpvalue<<std::endl;
        flag = false;
        
        std::cout<<"something wrong! decompress failed"<<std::endl;
        
      }

    }
    
    if(!flag){
        break;
    }
    */
    
  }
std::cout << "total compression rate:" << std::setprecision(4)
       << compressrate/blocks << std::endl;
std::cout<<"decompress all!"<<std::endl;
std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;
/*
std::cout<<"random access decompress!"<<std::endl;
std::cout << "random decoding time per int: " << std::setprecision(8)
     << randomaccesstime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "random decoding speed: " << std::setprecision(10)
     << data.size()/(randomaccesstime*1000) <<  std::endl;

    
*/
  
  //********************************************************************************
  // If you need to use differential coding, you can use
  // calls like these to get the deltas and recover the original
  // data from the deltas:

  
}
