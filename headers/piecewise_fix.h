/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */
#ifndef PIECEWISEFIX_H_
#define PIECEWISEFIX_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {




    
class piecewise_fix : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
  std::vector<lr> lrvec;
  std::vector<uint8_t> maxerror;


  




    
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *out, size_t nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];

    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    int *delta = new int[length];

    //lr.train(indexes, keys, length, 0.0001, 500);
    
    lr mylr;
    mylr.caltheta(indexes,keys,length);
    lrvec.push_back(mylr);
    
    //std::cout<<"Theta: "<<mylr.theta0<<" "<<mylr.theta1<<std::endl;
    
    free(indexes);
    free(keys);
    int max_error =0;
    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - (long long)(mylr.theta0+mylr.theta1*(double)i);
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
    /*
    if(in[0]==0){
    for(int i=0;i<length;i++){
        std::cout<<i<<" "<<delta[i]<<" "<<in[i]<<std::endl;
    }
    }
    */
    //std::cout<<"bit_length: "<<tmp_bit<<std::endl;
    if(tmp_bit>=31){
     out = write_delta_default(in,out,32,length);
    }
    else{
    out = write_delta(delta, out, tmp_bit, length);
    }
    maxerror.push_back(tmp_bit);
    
    free(delta);
    
    return out;
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
    //std::cout<<"decompressing all!"<<std::endl;
    lr mylr= lrvec[nvalue];
    if(maxerror[nvalue]>=31){
        read_all_default(in ,0,0, length, maxerror[nvalue],mylr.theta1,mylr.theta0, out);
    }
    else{
        read_all_bit_fix(in ,0,0, length, maxerror[nvalue],mylr.theta1,mylr.theta0, out);
    }

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    //double start =getNow();
    lr mylr= lrvec[nvalue];
    //double end = getNow();
    
    //double start2 = getNow();
    
    uint32_t tmp=0;
    if(maxerror[nvalue]>=31){
       tmp = read_bit_default(in ,maxerror[nvalue] , l, mylr.theta1,mylr.theta0, 0);
     }
    else{
       tmp = read_bit_fix(in ,maxerror[nvalue] , l, mylr.theta1,mylr.theta0, 0);
    }
    //double end2 = getNow();
    //std::cout<<"find lower bound time: "<<(end-start)*1000000000<<" read bit time: "<<(end2-start2)*1000000000<<" rate is: "<<(end-start)/(end2-start2)<<std::endl;
    return tmp;
    
   

}

uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *out,
                   size_t nvalue) {
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return out;
}
uint32_t *decodeArray( uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue) {
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return out;
}
uint32_t randomdecodeArray(uint32_t *in, const size_t l,uint32_t *out, size_t nvalue){
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return 1;
}
uint32_t get_block_nums(){
      return 1;
}    
std::string name() const {
    return "piecewise_fix"; 
}    
  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
