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
  int *method_times=new int[4];
  for(int i=0;i<4;i++){
            method_times[i]=0;
  }
  int totallines=0;
  std::ofstream outfile("../train_data/training.csv", std::ios::app);
  std::vector<std::string> codec_name={"piecewise_fix","nonlinear_fix","FOR","rle"};
  std::vector<std::string> data_name={
  "../train_data/training_piecewise.txt","../train_data/training_nonlinear.txt","../train_data/training_FOR.txt","../train_data/training_rle.txt"};
  std::vector<int> data_size ={393194,48122,64,1922};
  for (int pick_data=0;pick_data<data_name.size();pick_data++){
      std::ifstream srcFile(data_name[pick_data],std::ios::in); 
 
      if(!srcFile) { 
          std::cout << "error opening source file." << std::endl;
          return 0;
      }
      
      int len;
      uint32_t min;
      int logdelta;
      double quarter;
      double half;
      double threequarter;
      int num_distinct;
      int rl;
      int method;
      int line=0;
      std::string tmp;
      getline(srcFile,tmp);
      int sample_size = ceil(data_size[pick_data]/500);
      while(1){
          srcFile >> logdelta>>quarter>>half>>threequarter>>rl>>method;
          
          if(srcFile.eof()){break;}
          if(data_size[pick_data]<500){
                    outfile<< logdelta<<","<<quarter<<","<<half<<","<<threequarter<<","<<rl<<","<<method<<std::endl;
                    line++;
                    method_times[method]++;
          }
          else{
              if(random(sample_size)==0){
                    outfile<< logdelta<<","<<quarter<<","<<half<<","<<threequarter<<","<<rl<<","<<method<<std::endl;
                    line++;
                    method_times[method]++;
                }

          }
          


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
