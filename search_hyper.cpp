// searching for hyper-parameter like block number

#include "../headers/common.h"
#include "../headers/codecfactory.h"
#include "../headers/caltime.h"
#include "../headers/lr.h"
#include "../headers/search_hyper.h"
#include "../headers/create_feature.h"
#include "headers/microunit.h"
#include "headers/easylogging++.h"
#include "headers/MLP.h"
#include "headers/decision_tree.h"
#include "headers/file_manage.h"
#include "headers/model_selection.h"
using namespace Eigen;

INITIALIZE_EASYLOGGINGPP

using namespace Codecset;
int main()
{
    std::vector<std::string> codec_name = {"piecewise_fix","nonlinear_fix","FOR","rle"};
    std::vector<uint32_t> data;
    std::ifstream srcFile("../data/wf/newman.txt", std::ios::in);
    if (!srcFile)
    {
        std::cout << "error opening source file." << std::endl;
        return 0;
    }
    while (1)
    {

        uint32_t next;
        srcFile >> next;
        if (srcFile.eof())
        {
            break;
        }
        data.push_back(next);
    }
    srcFile.close();
    int N = data.size();

    if (data.size() == 0)
    {
        std::cout << "Empty vector" << std::endl;
        return 0;
    }
    std::cout << "vector size = " << data.size() << std::endl;
    std::cout << "vector size = " << data.size() * sizeof(uint32_t) / 1024.0 << "KB"
              << std::endl;
    int bsize[5] = {100, 200, 1000, 2000, 10000};
    int vote[5] ={0};
    int sample_time =1;
    double sample_rate = 0.001;
    int delta=0;
    double start = getNow();
    for(int i=0;i<(int)codec_name.size();i++){
        codec_vote tmp = pick_block_size(bsize,5,sample_time,sample_rate,N,data.data(),codec_name[i]);
        vote[tmp.select]+= 1.0/tmp.compression_rate;
    }
    int pos =0;
    int max_vote=0;
    for(int i=0;i<5;i++){
        if(vote[i]>max_vote){
            max_vote=vote[i];
            pos =i;
        }
    }
    int block_size = bsize[pos];
    
    double end = getNow();
    double search_time = end - start;

    
    int blocks = data.size() / block_size;
    if (blocks * block_size < N)
    {
        blocks++;
    } //handle with the last block, maybe < block_size
    std::cout << "Total blocks " << blocks << " block size " << block_size << std::endl;
    const int input_size = 5;
    const int number_classes = 4;

    const std::string weights = "../model_5dim.txt";
    std::ifstream infile(weights, std::ios::in);
    // prepare classifier
    DecisionTree model;
    model.rebuild(infile,number_classes);

//************************************ predict cocec ************************
  std::vector<IntegerCODEC*> codec_fac;
  for(int i=0;i<(int)codec_name.size();i++){
      IntegerCODEC &codec = *CODECFactory::getFromName(codec_name[i]);
      codec.init(blocks,block_size,delta);
      codec_fac.push_back(&codec);
  }
    
  std::vector<uint8_t*> block_start_vec;
  std::vector<int> method_vec;
  int totalsize = 0;
  
   start = getNow();
  double totaltime_realcom=0;
  for(int i=0;i<blocks;i++){

    seg_feature seg;
    seg.cal_feature(data.data()+(i*block_size),block_size);
    
    Eigen::MatrixXd tmp_feature = Eigen::MatrixXd::Zero(1 , input_size);
    tmp_feature<<seg.logdelta,seg.quarter,seg.half,seg.threequarter,seg.rl;
    Eigen::VectorXd pred(tmp_feature.rows());
    pred = model.predict( tmp_feature);
    size_t class_id=pred[0];
    
   
    double start2 = getNow();
    uint8_t * descriptor = (uint8_t*)malloc(block_size * sizeof(uint64_t)*2);
    uint8_t * res = descriptor;
    res = codec_fac[class_id]->encodeArray8(data.data()+(i*block_size),block_size ,descriptor,i);
    int tmp_size = (res-descriptor);
    double end2 = getNow();
    totaltime_realcom +=(end2-start2);
   
    method_vec.push_back(class_id);
    block_start_vec.push_back(descriptor);
    totalsize +=tmp_size;   
 
  }
    end = getNow();
   double totaltime = end -start;
   std::cout << "search time / compress time: " << std::setprecision(10) << search_time/totaltime<<  std::endl;
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

  bool flag = true;
  std::vector<uint32_t> recover(data.size());
  totaltime =0.0;
  std::cout<<"decompress all!"<<std::endl;
  start = getNow();
  for(int i=0;i<blocks;i++){
      //std::cout<<"block "<<(int)i<<" method "<<codec_name[method_vec[(int)i]]<<std::endl;
      codec_fac[method_vec[i]]->decodeArray8(block_start_vec[i], block_size, recover.data()+i*block_size, i);
      
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
  
      uint32_t tmpvalue = codec_fac[method_vec[(int)i/block_size]]->randomdecodeArray8(block_start_vec[(int)i/block_size], i%block_size, buffer.data(), i/block_size);
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
