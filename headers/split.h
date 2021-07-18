
#ifndef SPLIT_H_
#define SPLIT_H_

#include "common.h"
#include "Utils.h"
#include<math.h>
struct seg{//theta0+theta1*x
    int start;
    int end;
    double theta0;
    double theta1;
    int totalbyte;
    int tmpbit ;
    int* delta;
    
void caltheta(double x[], double y[], int m){
    double * real_x = new double[m];
    for(int i=0;i<m;i++){
        real_x[i]=i;
    }
    double avx= Utils::array_sum(real_x,m)*((double)1/m);
    double avy= Utils::array_sum(y,m)*((double)1/m);
    theta1=(Utils::array_sum(Utils::array_multiplication(real_x,y,m),m)-(double(m)*avx*avy))/(Utils::array_sum(Utils::array_multiplication(real_x,real_x,m),m)-(double(m)*avx*avx));
    theta0=avy - theta1*avx;
    
}
    
void calbit(double x[], double y[], int m){
    uint32_t max_error=0;
    delta = new int[m+20];
    for(int i=0;i<(long long)m;i++){
        int tmp = (long long) y[i] - (long long)(theta0 + (theta1*(double)(x[i]-x[0])));
        delta[i]=tmp;
        if(abs(tmp)>max_error){
            max_error = abs(tmp);
        }
    }
    int tmp_bit = 0;
    if(max_error > 0.01){
        tmp_bit = ceil(log2(max_error))+2;
    }
    else{
        tmp_bit = 2;
    }
    if(tmp_bit>=32){
        tmp_bit =32;
    }
    tmpbit=tmp_bit;
    totalbyte = 4+4+8+8+1+ceil(((double)tmp_bit*(double)m)/8.);
    
}


};
    






#endif 
