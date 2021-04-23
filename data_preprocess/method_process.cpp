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
  std::ofstream outfile("../train_data/training_rle.txt", std::ios::app);
  outfile<< "logdelta"<<"    "<<"0.25"<<"    "<<"0.5"<<"    "<<"0.75"<<"    "<<"rl"<<"    label"<<std::endl;
  std::vector<std::string> codec_name={"piecewise_fix","nonlinear_fix","FOR","rle"};
  std::vector<std::string> data_name={
  "../train_data/linear.txt","../train_data/normal.txt","../train_data/lognormal.txt","../train_data/books.txt","../train_data/fb.txt","../train_data/wiki.txt","../train_data/newman.txt"};

  for (int pick_data=0;pick_data<data_name.size();pick_data++){
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
      std::string tmp;
      getline(srcFile,tmp);
      while(1){
          srcFile >> logdelta>>quarter>>half>>threequarter>>rl>>method;
          
          if(srcFile.eof()){break;}
          if(method==2){
                outfile<< logdelta<<"    "<<quarter<<"    "<<half<<"    "<<threequarter<<"    "<<rl<<"    "<<method<<std::endl;
                line++;
          }
          


      }
      std::cout<<data_name[pick_data]<<" lines: "<<line<<std::endl;
      totallines+=line;
      srcFile.close();
    }
  
  std::cout<<"total lines : "<<totallines<<std::endl;


  outfile.close();
  
 
    




  
}
