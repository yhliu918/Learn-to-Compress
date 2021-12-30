
#ifndef LR_STRING_H_
#define LR_STRING_H_
#include "common.h"
#include "Utils.h"
#include<math.h>

typedef boost::multiprecision::mpf_float long_float;
struct string_lr{
    
long_float theta0 = 0;
long_float theta1 = 0;
long_float delta = 0;
    
void caltheta(std::vector<long_int>& x, std::vector<long_int>& y, int m){

    long_int sumx = 0;
    long_int sumy = 0;
    long_int sumxy = 0;
    long_int sumxx = 0;
    for(int i=0;i<m;i++){
        sumx = sumx + x[i];
        sumy = sumy + y[i];
        sumxx = sumxx+x[i]*x[i];
        sumxy = sumxy+x[i]*y[i];
    }
    long_float avx = sumx.convert_to<long_float>()/m;
    long_float avy = sumy.convert_to<long_float>()/m;
    long_float avxy = sumxy.convert_to<long_float>()/m;
    long_float avxx = sumxx.convert_to<long_float>()/m;

    long_float ccc= avxy - avx * avy ;
    long_float xxx = avxx - avx * avx;
    theta1= ccc/ xxx;
    theta0 = avy - theta1 * avx;
    
}
void caltheta_LOO(double x[], double y[], int m){
    double sumx= Utils::array_sum(x,m);
    double sumy= Utils::array_sum(y,m);
    double xy = Utils::array_sum(Utils::array_multiplication(x,y,m),m);
    double xx = Utils::array_sum(Utils::array_multiplication(x,x,m),m);
    for(int i=0;i<m;i++){
        double tmpavx = (sumx - x[i])*((double)1/(m-1));
        double tmpavy = (sumy - y[i])*((double)1/(m-1));
        double tmpxy  = xy - x[i]*y[i];
        double tmpxx  = xx - x[i]*x[i];
        theta1=(tmpxy -(double(m-1)*tmpavx*tmpavy))/(tmpxx-(double(m-1)*tmpavx*tmpavx));
        theta0=tmpavy - theta1*tmpavx;
        std::cout<<"Theta1: "<<theta1<<"Theta0: "<<theta0<<std::endl;
    }
    
    
}
    
};
    






#endif 
