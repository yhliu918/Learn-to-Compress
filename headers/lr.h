
#ifndef LR_H_
#define LR_H_
#include "bignumber.h"
#include "common.h"
#include "Utils.h"
#include<math.h>
struct lr{//theta0+theta1*x
    double theta0;
    double theta1;
    int delta;
    
void caltheta(double x[], double y[], int m){
    double avx= Utils::array_sum(x,m)*((double)1/m);
    double avy= Utils::array_sum(y,m)*((double)1/m);
    theta1=(Utils::array_sum(Utils::array_multiplication(x,y,m),m)-(double(m)*avx*avy))/(Utils::array_sum(Utils::array_multiplication(x,x,m),m)-(double(m)*avx*avx));
    theta0=avy - theta1*avx;
    
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
bool agree(double x, double y){
    double tmp = abs((x*theta1+theta0)-y);
    if(tmp<=(double)delta*delta){
        return true;
    }
    else{
        return false;
    }
}
    
};
    

struct lr_int{
    
int theta0 = 0;
int theta1 = 0;
    
void caltheta(uint32_t x[], uint32_t y[], int m){

    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumxx = 0;
    for(int i=0;i<m;i++){
        sumx = sumx + x[i];
        sumy = sumy + y[i];
        sumxx = sumxx+x[i]*x[i];
        sumxy = sumxy+x[i]*y[i];
    }
    
    double ccc= sumxy * m - sumx * sumy;
    double xxx = sumxx * m - sumx * sumx;

    theta1 = ccc/xxx;
    double theta1_d = ccc/xxx;
    theta0 = (sumy - theta1_d * sumx)/(double)m;
    
}
    
    
};


template <typename T>
inline uint32_t bits_int_T(T v)
{
  uint32_t r(0);
  int length = sizeof(T) * 8;
  if (length > 255 && v >= ((T)1 << (uint8_t)255))
  {
    v >>= 255;
    //v = v/2;
    r += 256;
  }
  if (length > 127 && v >= ((T)1 << (uint8_t)127))
  {
    v >>= 128;
    r += 128;
  }
  if (length > 63 && v >= ((T)1 << (uint8_t)63))
  {
    v >>= 64;
    r += 64;
  }
  if (length > 31 && v >= ((T)1 << (uint8_t)31))
  {
    v >>= 32;
    r += 32;
  }
  if (length > 15 && v >= ((T)1 << (uint8_t)15))
  {
    v >>= 16;
    r += 16;
  }
  if (length > 7 && v >= ((T)1 << (uint8_t)7))
  {
    v >>= 8;
    r += 8;
  }
  if (length > 3 && v >= ((T)1 << (uint8_t)3))
  {
    v >>= 4;
    r += 4;
  }
  if (length > 1 && v >= ((T)1 << (uint8_t)1))
  {
    v >>= 2;
    r += 2;
  }
  if (v >= (T)1)
  {
    r += 1;
  }
  
  return r;
}
    






#endif 
