#ifndef CREAT_FEATURE_H
#define CREAT_FEATURE_H

#include <vector>




struct seg_feature{
    int len;
    uint32_t avg;
    uint32_t max;
    uint32_t min;
    int num_distinct;
    int rl;//avg run length
    
    void cal_feature(uint32_t arr[], int length){
        max =arr[length-1];
        min = arr[0];
        uint64_t sum =0;
        int count =0;
        uint32_t last = arr[0];
        num_distinct = 1;
        for(int i=0;i<length;i++){
            if(arr[i]!=last){
                num_distinct++;
                last = arr[i];
            }
            sum +=arr[i];
        }
        avg = sum/length;
        len = length;
        rl = len/num_distinct;

    }
    
    void write_feature(std::ofstream &ff,int method){
        ff<< len <<"    "<<avg<<"    "<<min<<"    "<<max<<"    "<<num_distinct<<"    "<<rl<<"    "<<method<<std::endl;
        
    }
};



#endif 
