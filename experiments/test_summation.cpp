
#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"






int main() {
  using namespace Codecset;

  // We pick a CODEC
  IntegerCODEC &codec = *CODECFactory::getFromName("piecewise_fix");


  std::vector<uint32_t> data;

  double start = getNow();
  std::ifstream srcFile("../data/standard/books_200M_uint32.txt",std::ios::in); 
  while(1){
      
      uint32_t next ;
      srcFile >> next;
      if(srcFile.eof()){break;}
      data.push_back(next);

  }
  srcFile.close();
  int N = data.size();
  uint64_t sum =0;
  for(int i=0;i<N;i++){
    sum+=data[i];
  }
  double end = getNow();

  std::cout << "all summation time per int: " << std::setprecision(8)
     << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;


  int blocks =1000000;
  int block_size = data.size()/blocks;
  int delta =32;
  codec.init(blocks,block_size,delta);
  std::vector<uint8_t*> block_start_vec;
  std::vector<int> start_index;
  int totalsize = 0;
  
  std::ofstream outfile("save_descriptor.txt", std::ios::out);
  for(int i=0;i<blocks;i++){
    int block_length = block_size;
    if(i==blocks-1){
      block_length = N - (blocks-1)*block_size;
    }
    uint8_t * descriptor = (uint8_t*)malloc(block_length * sizeof(uint64_t)*2);
    uint8_t * res = descriptor;
    res = codec.encodeArray8(data.data()+(i*block_size),block_length ,descriptor,i);
    descriptor = (uint8_t*)realloc(descriptor, (res-descriptor));
    uint32_t* save_descriptor = reinterpret_cast<uint32_t*> (descriptor);
    int descriptor_size = (res-descriptor);
    int tmp_size = descriptor_size/4;
    if(tmp_size*4<descriptor_size){tmp_size++;}
    outfile<<descriptor_size<<std::endl;
    for(int j=0;j<tmp_size;j++){
      outfile<<save_descriptor[j]<<std::endl;
    }
    totalsize += (res-descriptor);
 
  }
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;
  


start = getNow();
  std::ifstream desFile("save_descriptor.txt",std::ios::in); 
  
  while(1){
      
      uint32_t next ;
      desFile >> next;
      int tmpsize = next/4;
      if(tmpsize*4<(int)next){tmpsize++;}
      uint32_t * descriptor = (uint32_t*)malloc(tmpsize*sizeof(uint32_t));
      for(int i=0;i<tmpsize;i++){
        desFile>>descriptor[i];
      }
      uint8_t* real_descriptor = reinterpret_cast<uint8_t*>(descriptor);
      real_descriptor = (uint8_t*)realloc(real_descriptor, next);
      if(desFile.eof()){break;}
      block_start_vec.push_back(real_descriptor);

  }
  desFile.close();
  std::cout<<"descriptor loaded"<<std::endl;

  std::cout<<"summation all!"<<std::endl;
  
  uint64_t sum2 =0;
  
  for(int i=0;i<blocks;i++){
      //sum2+= codec.summation(block_start_vec[i],block_size,block_size);
  }
  end = getNow();

  if(sum!=sum2){
    std::cout<<"summation failed!"<<std::endl;
  }
  

std::cout << "all summation time per int: " << std::setprecision(8)
     << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;

  
   for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }



  
}
