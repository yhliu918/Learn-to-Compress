
#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"






int main() {
  using namespace Codecset;

  // We pick a CODEC
  IntegerCODEC &codec = *CODECFactory::getFromName("delta_my");


  std::vector<uint32_t> data;
  std::ifstream srcFile("../data/fb/fb-289000.txt",std::ios::in); 
  while(1){
      
      uint32_t next ;
      srcFile >> next;
      if(srcFile.eof()){break;}
      data.push_back(next);

  }
  srcFile.close();
  int N = data.size();

  std::vector<bool> bitmap;
  std::ifstream bitFile("../data/bitmap/bitmap_cluster_0.5_1000_289000.txt",std::ios::in); 
  while(1){
      
      bool next ;
      bitFile >> next;
      if(bitFile.eof()){break;}
      bitmap.push_back(next);

  }
  bitFile.close();

  int blocks =1000;
  int block_size = data.size()/blocks;
  blocks = data.size()/block_size;
  if(blocks*block_size<N){blocks++;} //handle with the last block, maybe < block_size
  std::cout<<"Total blocks "<<blocks<<" block size "<<block_size<<std::endl;
  int delta =32;
  codec.init(blocks,block_size,delta);
  std::vector<uint8_t*> block_start_vec;
  std::vector<int> start_index;
  int totalsize = 0;
  for(int i=0;i<blocks;i++){
    int block_length = block_size;
    if(i==blocks-1){
      block_length = N - (blocks-1)*block_size;
    }
    uint8_t * descriptor = (uint8_t*)malloc(block_length * sizeof(uint64_t));
    uint8_t * res = descriptor;
    res = codec.encodeArray8(data.data()+(i*block_size),block_length ,descriptor,i);
    descriptor = (uint8_t*)realloc(descriptor, (res-descriptor));
    block_start_vec.push_back(descriptor);
    totalsize += (res-descriptor);
 
  }
  
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;

/*
  uint32_t* buffer = NULL;
  double start = getNow();
  uint32_t tmpvalue =0;
  for(int i=0;i<N;i++){    
    if(bitmap[i]){
      tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)i/block_size], i%block_size, buffer, i/block_size);
    }
  }
  double end = getNow();
  std::cout << "access time per int: " << std::setprecision(8)
     << (end - start)/data.size()  * 1000000000 << "ns" << std::endl;

*/

/*
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

   for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }



  
}
