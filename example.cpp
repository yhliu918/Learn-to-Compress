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
#include "../headers/create_feature.h"



int main() {
  using namespace Codecset;

  // We pick a CODEC


  std::vector<uint32_t> data;
  std::ifstream srcFile("../data/wf/newman.txt",std::ios::in); 
  std::ofstream outfile("out.txt", std::ios::app);
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
 
    
  int blocks =100;
  int block_size = data.size()/blocks;
  std::vector<std::string> codec_name={"piecewise_fix","nonlinear_fix","FOR"};
  int delta =0;
  
  std::vector<uint8_t*> block_start_vec;
  std::vector<int> method_vec;
  int totalsize = 0;
  //outfile<< "len" <<"    "<<"avg"<<"    "<<"min"<<"    "<<"max"<<"    "<<"num_distinct"<<"    "<<"rl"<<"    label"<<std::endl;
  for(int i=0;i<blocks;i++){
    int min_size = block_size * 8;
    int method =0;
    uint8_t * tmp_des = (uint8_t*)malloc(block_size * sizeof(uint64_t)*2);
    seg_feature seg;
    seg.cal_feature(data.data()+(i*block_size),block_size);
    
    for(int j=0;j<(int)codec_name.size();j++){
      IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[j]);
      codec.init(blocks,block_size,delta);
      
      uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t)*2);
      uint8_t * res = descriptor;
      res = codec.encodeArray8(data.data()+(i*block_size),block_size ,descriptor,i);
      int tmp_size = (res-descriptor);
      if(tmp_size<min_size){
              min_size = tmp_size;
              method = j;
              memcpy(tmp_des,descriptor,tmp_size);
              tmp_des = (uint8_t*)realloc(tmp_des, tmp_size);
              free(descriptor);
      }
      
    }
   //seg.write_feature(outfile,method);
   method_vec.push_back(method);
   block_start_vec.push_back(tmp_des);
   totalsize +=min_size;   
 
  }
  outfile.close();

  /*
  for(int i=0;i<blocks;i++){
      std::cout<<"block "<<i<<" method "<<codec_name[method_vec[i]]<<std::endl;
  }
  */
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;

  bool flag =true;
  std::vector<uint32_t> recover(data.size());
  double totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
   double start = getNow();
  for(int i=0;i<blocks;i++){
      //std::cout<<"block "<<(int)i<<" method "<<codec_name[method_vec[(int)i]]<<std::endl;
      IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[method_vec[i]]);
      codec.init(blocks,block_size,delta);
      codec.decodeArray8(block_start_vec[i], block_size, recover.data()+i*block_size, i);
      
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

  std::cout<<"random access decompress!"<<std::endl; 
  std::vector<uint32_t> buffer(data.size());
  double randomaccesstime =0.0;
   start = getNow();
   uint32_t mark=0;
    
  for(int i=0;i<N;i++){ 
      IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[method_vec[(int)i/block_size]]);
      codec.init(blocks,block_size,delta);
      //std::cout<<"block "<<(int)i/block_size<<" method "<<codec_name[method_vec[(int)i/block_size]]<<std::endl;
      uint32_t tmpvalue = codec.randomdecodeArray8(block_start_vec[(int)i/block_size], i%block_size, buffer.data(), i/block_size);
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
    
   for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }



  
}
