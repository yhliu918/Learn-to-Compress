
#ifndef LR_H_
#define LR_H_
#include <math.h>
#include "Utils.h"
#include "bignumber.h"
#include "common.h"
#include "string/leco_uint256.h"
struct lr {  // theta0+theta1*x
  double theta0;
  double theta1;
  int delta;


  void caltheta(std::vector<double>& x, std::vector<double>& y, int m) {
    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumxx = 0;
    for (int i = 0; i < m; i++) {
      sumx = sumx + x[i];
      sumy = sumy + y[i];
      sumxx = sumxx + x[i] * x[i];
      sumxy = sumxy + x[i] * y[i];
    }

    double ccc = sumxy * m - sumx * sumy;
    double xxx = sumxx * m - sumx * sumx;

    theta1 = ccc / xxx;
    theta0 = (sumy - theta1 * sumx) / (double)m;
  }

  void caltheta_LOO(double x[], double y[], int m) {
    double sumx = Utils::array_sum(x, m);
    double sumy = Utils::array_sum(y, m);
    double xy = Utils::array_sum(Utils::array_multiplication(x, y, m), m);
    double xx = Utils::array_sum(Utils::array_multiplication(x, x, m), m);
    for (int i = 0; i < m; i++) {
      double tmpavx = (sumx - x[i]) * ((double)1 / (m - 1));
      double tmpavy = (sumy - y[i]) * ((double)1 / (m - 1));
      double tmpxy = xy - x[i] * y[i];
      double tmpxx = xx - x[i] * x[i];
      theta1 = (tmpxy - (double(m - 1) * tmpavx * tmpavy)) /
               (tmpxx - (double(m - 1) * tmpavx * tmpavx));
      theta0 = tmpavy - theta1 * tmpavx;
      std::cout << "Theta1: " << theta1 << "Theta0: " << theta0 << std::endl;
    }
  }
  bool agree(double x, double y) {
    double tmp = abs((x * theta1 + theta0) - y);
    if (tmp <= (double)delta * delta) {
      return true;
    } else {
      return false;
    }
  }
};

template <typename T>
struct lr_int_T{//theta0+theta1*x
    double theta0;
    double theta1;

    
void caltheta(const T *y, int m){

    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumxx = 0;
    for(int i=0;i<m;i++){
        sumx = sumx + (double)i;
        sumy = sumy + (double)y[i];
        sumxx = sumxx+(double)i*i;
        sumxy = sumxy+(double)i*y[i];
    }
    
    double ccc= sumxy * m - sumx * sumy;
    double xxx = sumxx * m - sumx * sumx;

    theta1 = ccc/xxx;
    theta0 = (sumy - theta1 * sumx)/(double)m;
    
}

};

struct lr_int{//theta0+theta1*x
    double theta0;
    double theta1;
  
  void caltheta(uint32_t *y, int m){

    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumxx = 0;
    for(int i=0;i<m;i++){
        sumx = sumx + (double)i;
        sumy = sumy + (double)y[i];
        sumxx = sumxx+(double)i*i;
        sumxy = sumxy+(double)i*y[i];
    }
    
    double ccc= sumxy * m - sumx * sumy;
    double xxx = sumxx * m - sumx * sumx;

    theta1 = ccc/xxx;
    theta0 = (sumy - theta1 * sumx)/(double)m;
    
}

};

// assume v > 0
template <typename T,
          typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
inline uint32_t bits_int_T(T v) {
  if(v<0){
    v=-v;
  }
  if (v == 0) return 0;
#if defined(__clang__) || defined(__GNUC__)
  // std::cout<<__builtin_clzll(v)<<" "<<64 - __builtin_clzll(v)<<std::endl;
  return 64 - __builtin_clzll(v);
#else
  assert(false);
#endif
}
template <typename T, typename std::enable_if<std::is_same<__uint128_t, T>::value ||
                                                  std::is_same<leco_uint256, T>::value,
                                              bool>::type = true>
inline uint32_t bits_int_T(T v) {
  if(v<0){
    v = -v;
  }
  uint32_t r(0);
  constexpr int length = sizeof(T) * 8;
  if constexpr (length > 255) {
    if (v >= ((T)1 << (uint8_t)255)) {
      v >>= 255;
      // v = v/2;
      r += 256;
    }
  }
  if constexpr (length > 127) {
    if (v >= ((T)1 << (uint8_t)127)) {
      v >>= 128;
      r += 128;
    }
  }
  if constexpr (length > 63) {
    if (v >= ((T)1 << (uint8_t)63)) {
      v >>= 64;
      r += 64;
    }
  }
  if (length > 31 && v >= ((T)1 << (uint8_t)31)) {
    v >>= 32;
    r += 32;
  }
  if (length > 15 && v >= ((T)1 << (uint8_t)15)) {
    v >>= 16;
    r += 16;
  }
  if (length > 7 && v >= ((T)1 << (uint8_t)7)) {
    v >>= 8;
    r += 8;
  }
  if (length > 3 && v >= ((T)1 << (uint8_t)3)) {
    v >>= 4;
    r += 4;
  }
  if (length > 1 && v >= ((T)1 << (uint8_t)1)) {
    v >>= 2;
    r += 2;
  }
  if (v >= (T)1) {
    r += 1;
  }

  return r;
}

#endif
