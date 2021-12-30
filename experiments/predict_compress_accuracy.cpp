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
  std::ifstream srcFile("../data/standard/lognormal_200M_uint32.txt",std::ios::in); 
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
  
	
  
  int block_size =20000;
  int blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } //handle with the last block, maybe < block_size
  int delta =0;
  
  std::vector<IntegerCODEC*> codec_fac;
  std::vector<std::string> codec_name={"piecewise_fix","FOR","rle"};
  //std::vector<std::string> codec_name={"piecewise_fix"};
  for(int i=0;i<(int)codec_name.size();i++){
      IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[i]);
      codec.init(blocks,block_size,delta);
      codec_fac.push_back(&codec);
  }
    
  std::vector<int> method_vec;
  int totalsize = 0;
  //outfile<< "len" <<"    "<<"avg"<<"    "<<"min"<<"    "<<"max"<<"    "<<"num_distinct"<<"    "<<"rl"<<"    label"<<std::endl;
  double start = getNow();
  double totaltime_realcom=0;
  double percent = 1/(double)blocks;
//******************* PREDICT **************************************
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
    uint8_t * descriptor = (uint8_t*)malloc(block_length* sizeof(uint64_t)*2);
    uint8_t * res = descriptor;
    res = codec_fac[pick_method]->encodeArray8(data.data()+(i*block_size),block_length ,descriptor,i);
    int tmp_size = (res-descriptor);
    free(descriptor);
   //seg.write_feature(outfile,method);
    method_vec.push_back(pick_method);
    totalsize +=tmp_size;   
 
  }
  //outfile.close();
   double end = getNow();
   double compressrate = (totalsize)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate << std::endl;

  //******************* EXHAUSTIVE **************************************
  std::vector<int> method_vec_truth;
  int totalsize_best = 0;
  for(int i=0;i<blocks;i++){
    int min_size = block_size * 8;
    int method =0;
    uint8_t * tmp_des = (uint8_t*)malloc(block_size * sizeof(uint64_t)*2);
    //seg_feature seg;
    //seg.cal_feature(data.data()+(i*block_size),block_size);
    
    for(int j=0;j<(int)codec_name.size();j++){
      int block_length = block_size;
      if(i==blocks-1){
        block_length = N - (blocks-1)*block_size;
      }
      uint8_t * descriptor = (uint8_t*)malloc(block_length * sizeof(uint64_t)*2);
      uint8_t * res = descriptor;
      res = codec_fac[j]->encodeArray8(data.data()+(i*block_size),block_length ,descriptor,i);
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
   method_vec_truth.push_back(method);
   totalsize_best +=min_size;   
 
  }
  double compressrate_best = (totalsize_best)*100.0  / (4*N*1.0);
  std::cout << "total compression rate:" << std::setprecision(4)<< compressrate_best << std::endl;

  int *times_correct= new int[codec_name.size()];
  int *times= new int[codec_name.size()];
  for(int i=0;i<(int)codec_name.size();i++){
      times[i]=0;
      times_correct[i]=0;
  }
  int totalcorrect = 0;
  for(int i=0;i<blocks;i++){
      times[method_vec_truth[i]]++;
      if(method_vec_truth[i] == method_vec[i]){
        times_correct[method_vec_truth[i]]++;
        totalcorrect++;
      }

  }
  for(int i=0;i<(int)codec_name.size();i++){
      std::cout<< "method "<<codec_name[i]<<" correct percentage "<<(double)times_correct[i]/(double)times[i]<<std::endl;
  }
  std::cout<< "total correct "<<(double)totalcorrect/(double)blocks<<std::endl;
}
