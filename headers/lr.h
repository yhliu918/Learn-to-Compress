/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */
#ifndef LR_H_
#define LR_H_

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
    






#endif 
