
#ifndef LR_STRING_H_
#define LR_STRING_H_
#include "common.h"
#include "Utils.h"
#include<math.h>

typedef boost::multiprecision::mpf_float long_float;
struct string_lr{
    
long_int theta0 = 0;
long_int theta1 = 0;
    
void caltheta(std::vector<int>& x, std::vector<long_int>& y, int m){

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
    
    long_int ccc= sumxy * m - sumx * sumy;
    long_int xxx = sumxx * m - sumx * sumx;

    theta1 = ccc/xxx;
    theta0 = (sumy - theta1 * sumx)/m;
    
}
    
    
};
    



#endif 
