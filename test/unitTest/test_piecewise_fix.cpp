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
  IntegerCODEC &codec = *CODECFactory::getFromName("piecewise_fix");
  
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

  size_t N = 200000000;
  int blocks =100000;
  std::vector<uint32_t> data(N);
  FILE *fpRead=fopen("../data/linear_200M_uint32.txt","r");
  
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
  bool flag = true;
  std::vector<int> start_index;
  std::vector<uint8_t> compdata(4*data.size() + 1024); 
  uint8_t* out = compdata.data();
  uint8_t* ind = compdata.data();
  
  for(int i=0;i<blocks;i++){

    out = codec.encodeArray8(data.data()+(i*block_size),block_size, out ,N);
    start_index.push_back(ind - compdata.data());
    ind = out;
    
  }
  int totalsize = ind - compdata.data()+ blocks*(4+8+8+1);
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;
    
  std::vector<uint32_t> recover(data.size());
  double totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
  for(int i=0;i<blocks;i++){
      double start = getNow();
      codec.decodeArray8(compdata.data()+start_index[i], block_size, recover.data()+i*block_size, i);
      double end = getNow();
      totaltime += (end - start);
      for(int j=0;j<block_size;j++){
        if(data[j+i*block_size]!=recover[j+i*block_size]){
          std::cout<<"block: "<<i<<"num: "<<j<< "true is: "<<data[j+i*block_size]<<" predict is: "<<recover[j+i*block_size]<<std::endl;
          std::cout<<"something wrong! decompress failed"<<std::endl;
          flag = false;
          break;
         }
         
       }
       if(!flag){
          break;
       }
  }


   
  std::cout<<"random access decompress!"<<std::endl; 
  std::vector<uint32_t> buffer(data.size());
  double randomaccesstime =0.0;
  for(int i=0;i<N;i++){    
      double start = getNow();
      uint32_t tmpvalue = codec.randomdecodeArray8(compdata.data()+start_index[i/block_size], i%block_size, buffer.data(), i/block_size);
      double end = getNow();
      randomaccesstime+=(end-start);

      //std::cout<<"processing...  "<<j<<" / "<<N<<std::endl;

      if(data[i]!=tmpvalue){
        
        std::cout<<"num: "<<i<< "true is: "<<data[i]<<" predict is: "<<tmpvalue<<std::endl;
        flag = false;
        std::cout<<"something wrong! decompress failed"<<std::endl;
        
      }
    if(!flag){
        break;
    }
      
  }



std::cout<<"decompress all!"<<std::endl;
std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;
    
std::cout<<"random access decompress!"<<std::endl;
std::cout << "random decoding time per int: " << std::setprecision(8)
     << randomaccesstime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "random decoding speed: " << std::setprecision(10)
     << data.size()/(randomaccesstime*1000) <<  std::endl;



  
  //********************************************************************************
  // If you need to use differential coding, you can use
  // calls like these to get the deltas and recover the original
  // data from the deltas:

  
}
