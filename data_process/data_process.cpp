// processing training data

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
  int *method_times=new int[4];
  for(int i=0;i<4;i++){
            method_times[i]=0;
  }
  int totallines=0;
  std::ofstream outfile("../train_data_try/training_largesample.csv", std::ios::out);
  std::vector<std::string> codec_name={"piecewise_fix","FOR","rle"};
  std::vector<std::string> data_name={
  "../train_data_try/training_piecewise.txt","../train_data_try/training_FOR.txt","../train_data_try/training_rle.txt"};
  std::vector<int> data_size ={1332195,1332200,1332200};
  for (int pick_data=0;pick_data<(int)data_name.size();pick_data++){
      std::ifstream srcFile(data_name[pick_data],std::ios::in); 
 
      if(!srcFile) { 
          std::cout << "error opening source file." << std::endl;
          return 0;
      }
      

      int logdelta;
      double quarter;
      double half;
      double threequarter;
      int rl;
      int method;
      int line=0;
      double percent;
      double compressrate;
      std::string tmp;
      getline(srcFile,tmp);
      int sample_size = ceil(data_size[pick_data]/2000);
      while(1){
          srcFile >> logdelta>>quarter>>half>>threequarter>>rl>>method>>percent>>compressrate;
          //outfile<<logdelta<<","<<quarter<<","<<half<<","<<threequarter<<","<<rl<<","<<method<<","<<percent<<","<<compressrate<<std::endl;
          if(srcFile.eof()){break;}
          
          if(data_size[pick_data]<500){
                    outfile<<logdelta<<","<<quarter<<","<<half<<","<<threequarter<<","<<rl<<","<<method<<","<<percent<<","<<compressrate<<std::endl;
                    line++;
                    method_times[method]++;
          }
          else{
              if(random(sample_size)==0){
                    outfile<<logdelta<<","<<quarter<<","<<half<<","<<threequarter<<","<<rl<<","<<method<<","<<percent<<","<<compressrate<<std::endl;
                    line++;
                    method_times[method]++;
                }

          }
          
        //line++;
        //method_times[method]++;
          


      }
      std::cout<<data_name[pick_data]<<" lines: "<<line<<std::endl;
      totallines+=line;
      srcFile.close();
    }
  


  for(int i=0;i<4;i++){
      std::cout<<codec_name[i]<<" "<<"percent "<<(double)method_times[i]/totallines<<std::endl;

  }

  outfile.close();
  
 
    




  
}
