#ifndef CREAT_FEATURE_H
#define CREAT_FEATURE_H

#include <vector>




struct seg_feature{
    int logdelta;
    double quarter;
    double half;
    double threequarter;
    int rl;//avg run length
    
    void cal_feature(uint32_t arr[], int length){
        double max =arr[length-1];
        double min = arr[0];
        double delta = max-min;
        if(delta>0.01){
            logdelta=ceil(log2(delta))+2;
        }
        else{
            logdelta=2;
        }
        uint32_t last = arr[0];
        int num_distinct = 1;
        for(int i=0;i<length;i++){
            if(i==length/4){
                quarter = (double)(arr[i]-min)/(double)(max-min+0.0001);
                quarter = fabs((double)(quarter*10.0-2.5)*10.0);
            }
            if(i==length/2){
                half = (double)(arr[i]-min)/(double)(max-min+0.0001);
                half = fabs((double)(half*10.0-5.0)*10.0);
            }
            if(i==length*3/4){
                threequarter = (double)(arr[i]-min)/(double)(max-min+0.0001);
                threequarter= fabs((double)(threequarter*10.0-7.5)*10.0);
            }
            if(arr[i]!=last){
                num_distinct++;
                last = arr[i];
            }
        }
        int len = length;
        rl = len/num_distinct;
        

    }
    
    void write_feature(std::ofstream &ff,int method, double percent, double compressrate){
        ff<<logdelta<<"    "<<quarter<<"    "<<half<<"    "<<threequarter<<"    "<<rl<<"    "<<method<<"    "<<percent<<"    "<<compressrate<<std::endl;
        
    }
};



#endif 
