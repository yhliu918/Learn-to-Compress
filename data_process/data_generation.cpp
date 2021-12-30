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

int random(int m){
        return rand()%m;
}


int main() {
  using namespace Codecset;

  // We pick a CODEC
  std::ofstream outfile("../regression_model_train_data/newman.txt", std::ios::out);
  //outfile<< "logdelta"<<"    "<<"0.25"<<"    "<<"0.5"<<"    "<<"0.75"<<"    "<<"outlie"<<"    "<<"num_distinct"<<"    "<<"rl"<<"    std_delta"<<"    label"<<std::endl;
  std::vector<std::string> codec_name={"piecewise_fix","FOR","rle"};
  std::vector<std::string> data_name={"/home/ssq/Learn-to-Compress/data/linear_200M_uint32.txt","/home/ssq/Learn-to-Compress/data/normal_200M_uint32.txt","/home/ssq/Learn-to-Compress/data/lognormal_200M_uint32.txt","/home/ssq/Learn-to-Compress/data/books_200M_uint32.txt",
  "/root/Learn-to-Compress/data/gap/gap_47_0.8_100M_uint32.txt","/home/ssq/Learn-to-Compress/data/fb/fb-289000.txt","/home/ssq/Learn-to-Compress/data/wf/wiki.txt","/home/ssq/Learn-to-Compress/data/wf/newman.txt"};
  std::vector<int> max_block={100000,100000,100000,1000000,100000, 1000,1000,1000};
  std::vector<int> total_size={200000000,200000000,200000000,200000000,100000000, 289000,2076000,233000};
  for (int pick_data=6;pick_data<7;pick_data++){
      std::vector<uint32_t> data;
      std::ifstream srcFile(data_name[pick_data],std::ios::in); 
 
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
      int init_block = 100;
      if(pick_data<5){
          init_block=10000;
      }
      
      
      while(init_block<=max_block[pick_data]){
          int blocks =init_block;
          std::cout<<"data: "<<data_name[pick_data]<<" blocks "<< blocks<<std::endl;
          int block_size = data.size()/blocks;
          int delta =0;
    
          std::vector<IntegerCODEC*> codec_fac;
  
          for(int i=0;i<(int)codec_name.size();i++){
              IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[i]);
              codec.init(blocks,block_size,delta);
              codec_fac.push_back(&codec);
          }
    
          std::vector<uint8_t*> block_start_vec;
          std::vector<int> method_vec;
      

          for(int i=0;i<blocks;i++){
            int min_size = block_size * 8;
            int method =0;
            seg_feature seg;
            seg.cal_feature(data.data()+(i*block_size),block_size);
    
            for(int j=0;j<(int)codec_name.size();j++){
      
              uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t)*2);
              uint8_t * res = descriptor;
              res = codec_fac[j]->encodeArray8(data.data()+(i*block_size),block_size ,descriptor,i);
              double tmp_size = (res-descriptor)/(block_size*4.0);
              double percent = (double)block_size/(double)total_size[pick_data];
              seg.write_feature(outfile,j,percent,tmp_size);
              free(descriptor);
            }
            
          }
          
          init_block *=10;

      }
      

 
  }
  


  

  outfile.close();
  
 
    




  
}
