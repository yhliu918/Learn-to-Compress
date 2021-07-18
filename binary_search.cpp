// test binary search

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"

using namespace Codecset;
IntegerCODEC &codec = *CODECFactory::getFromName("FOR");
std::vector<uint32_t> data;
std::vector<uint8_t*> block_start_vec;
int blocks =1000000;
int block_size = 0;
uint32_t* buffer = (uint32_t*)malloc(200000000 * sizeof(uint32_t));
bool Bsearch( int low, int high, uint32_t key)
{
    int mid;
    if (low > high)
    {
        return false;
    }
    mid = (low + high) / 2;
    if (data[mid] == key)
    {
         return true;
     }
     else if (data[mid] > key)
     {
        return  Bsearch(low, mid - 1, key);
     }
     else if (data[mid] < key)
     {
        return  Bsearch(mid + 1, high, key);
     }
}

bool ourBsearch( int low, int high, uint32_t key)
{
    int mid;
    if (low > high)
    {
        return false;
    }
    mid = (low + high) / 2;
    uint32_t data_mid = codec.randomdecodeArray8(block_start_vec[(int)mid/block_size], mid%block_size, buffer, mid/block_size);
    if (data_mid == key)
    {
         return true;
     }
     else if (data_mid > key)
     {
        return  ourBsearch(low, mid - 1, key);
     }
     else if (data_mid < key)
     {
        return  ourBsearch(mid + 1, high, key);
     }
}


int main() {
  

  // We pick a CODEC
  std::ifstream srcFile("../data/books_200M_uint32.txt",std::ios::in); 
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
 
    
  
  int delta =32;
  block_size = data.size()/blocks;
  codec.init(blocks,block_size,delta);
  std::vector<int> start_index;
  int totalsize = 0;
  for(int i=0;i<blocks;i++){
    uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t));
    uint8_t * res = descriptor;
    res = codec.encodeArray8(data.data()+(i*block_size),block_size ,descriptor,i);
    descriptor = (uint8_t*)realloc(descriptor, (res-descriptor));
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

  }
      double end = getNow();
      totaltime += (end - start);

std::cout << "all decoding time per int: " << std::setprecision(8)
     << totaltime / data.size() * 1000000000 << "ns" << std::endl;
std::cout << "all decoding speed: " << std::setprecision(10)
     << data.size()/(totaltime*1000) <<  std::endl;



  int sample_size = 100000;
  start = getNow();
  for(int i=0;i<sample_size;i++){
      uint32_t tmpkey = rand();
      //std::cout<<tmpkey<<std::endl;
      ourBsearch(0,N,tmpkey);
      
  }
  end = getNow();
  double ourbinarytime = end - start;
  std::cout << "binary time per time: " << std::setprecision(8)
     << ourbinarytime / sample_size * 1000000000 << "ns" << std::endl;

   for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }



  
}
