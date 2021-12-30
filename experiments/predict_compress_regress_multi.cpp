// selecting codec by predict

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
#include "../headers/create_feature.h"
#include "../headers/microunit.h"
#include "../headers/easylogging++.h"
#include "../headers/MLP.h"
#include "../headers/regress_tree.h"
#include "../headers/file_manage.h"
#include "../headers/model_selection.h"
using namespace Eigen;

INITIALIZE_EASYLOGGINGPP

const int input_size = 7;

std::vector<std::string> weights = {"../reg_model/reg_model_piecewise.txt","../reg_model/reg_model_FOR.txt","../reg_model/reg_model_rle.txt"};


int main() {
  using namespace Codecset;

  // We pick a CODEC


  std::vector<uint32_t> data;
  std::ifstream srcFile("../data/standard/normal_200M_uint32.txt",std::ios::in); 
  //std::ofstream outfile("out.txt", std::ios::app);
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
 
  // prepare classifier
  
  std::vector<RegressionTree> models;
  for (int i=0;i<(int)weights.size();i++){
    std::ifstream infile(weights[i], std::ios::in);
    RegressionTree model;
    model.rebuild(infile,0);
    models.push_back(model);
    infile.close();
  }
  
	
  
  int blocks =1000;
  int block_size = data.size()/blocks;
  int delta =0;
  
  std::vector<IntegerCODEC*> codec_fac;
  std::vector<std::string> codec_name={"piecewise_fix","FOR","rle"};
  //std::vector<std::string> codec_name={"piecewise_fix"};
  for(int i=0;i<(int)codec_name.size();i++){
      IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[i]);
      codec.init(blocks,block_size,delta);
      codec_fac.push_back(&codec);
  }
    
  std::vector<uint8_t*> block_start_vec;
  std::vector<int> method_vec;
  int totalsize = 0;
  //outfile<< "len" <<"    "<<"avg"<<"    "<<"min"<<"    "<<"max"<<"    "<<"num_distinct"<<"    "<<"rl"<<"    label"<<std::endl;
  double start = getNow();
  double totaltime_realcom=0;
  double percent = 1/(double)blocks;
  for(int i=0;i<blocks;i++){
  int block_length = block_size;
    if(i==blocks-1){
      block_length = N - (blocks-1)*block_size;
    }
    seg_feature seg;
    seg.cal_feature(data.data()+(i*block_size),block_length);
    int pick_method =0;
    double pick_rate = 1.0;
    for(int j=0;j<(int)codec_name.size();j++){
        Eigen::MatrixXd tmp_feature = Eigen::MatrixXd::Zero(1 , input_size);
        tmp_feature<<seg.logdelta,seg.quarter,seg.half,seg.threequarter,seg.rl,j,percent;
        VectorXd pred(tmp_feature.rows());
        pred = models[j].predict( tmp_feature);
        double pred_rate = pred[0];
        //std::cout<<"method "<< codec_name[j]<<" pred rate "<<pred_rate<<std::endl;
        if(pred_rate<pick_rate){
          pick_rate = pred_rate;
          pick_method = j;
        }
    }
   
    double start2 = getNow();
    uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t)*2);
    uint8_t * res = descriptor;
    
    res = codec_fac[pick_method]->encodeArray8(data.data()+(i*block_size),block_length ,descriptor,i);
    int tmp_size = (res-descriptor);
    double end2 = getNow();
    totaltime_realcom +=(end2-start2);
   //seg.write_feature(outfile,method);
    method_vec.push_back(pick_method);
    block_start_vec.push_back(descriptor);
    totalsize +=tmp_size;   
 
  }
  //outfile.close();
   double end = getNow();
   double totaltime = end -start;
   std::cout << "compress speed: " << std::setprecision(10) << data.size()/(totaltime*1000) <<  std::endl;
   std::cout << "real compress speed: " << std::setprecision(10) << data.size()/(totaltime_realcom*1000) <<  std::endl;
  /*
  for(int i=0;i<blocks;i++){
      std::cout<<"block "<<i<<" method "<<codec_name[method_vec[i]]<<std::endl;
  }
  */
  int *times= new int[codec_name.size()];
  for(int i=0;i<(int)codec_name.size();i++){
      times[i]=0;
  }
  for(int i=0;i<blocks;i++){
      times[method_vec[i]]++;
  }
  for(int i=0;i<(int)codec_name.size();i++){
      std::cout<< "method "<<codec_name[i]<<" percentage "<<(double)times[i]/(double)blocks<<std::endl;
  }
  double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;

  std::vector<uint32_t> recover(data.size());
  totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
  start = getNow();
  for(int i=0;i<blocks;i++){
    int block_length = block_size;
      if(i==blocks-1){
        block_length = N - (blocks-1)*block_size;
      }
      //std::cout<<"block "<<(int)i<<" method "<<codec_name[method_vec[(int)i]]<<std::endl;
      codec_fac[method_vec[i]]->decodeArray8(block_start_vec[i], block_length, recover.data()+i*block_size, i);
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
   end = getNow();
  totaltime = (end - start);

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
      int block_length = block_size;
        if (i == blocks - 1)
        {
            block_length = N - (blocks - 1) * block_size;
        }
      uint32_t tmpvalue = codec_fac[method_vec[(int)i/block_size]]->randomdecodeArray8(block_start_vec[(int)i/block_size], i%block_size, buffer.data(), i/block_size);
       mark+=tmpvalue;
      
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
    
   for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }



  
}
