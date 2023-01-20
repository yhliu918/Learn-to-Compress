
#ifndef SPLINE_FIX_H_
#define SPLINE_FIX_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "spline_lr.h"
#define INF 0x7f7fffff

namespace Codecset {




    
class spline_fix : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
  using IntegerCODEC::init;


  int block_num;
  int block_size;
  
void init(int blocks, int blocksize,int extra){
      block_num=blocks;
      block_size=blocksize;
      
}



// bit + theta0 + theta1 + delta   
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];
    uint8_t *out = res;
    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    int *delta = new int[length];

    //lr.train(indexes, keys, length, 0.0001, 500);
    
    spline_lr mylr;
    mylr.caltheta(indexes,keys,length);
    
    //std::cout<<"Theta: "<<mylr.theta0<<" "<<mylr.theta1<<std::endl;
    
    free(indexes);
    free(keys);
    int max_error =0;
    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - (long long)(mylr.alpha+mylr.theta1*(double)i+mylr.theta2*(double)i*(double)i+mylr.theta3*(double)i*(double)i*(double)i);
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
    out[0]=(uint8_t)tmp_bit;
    out++;
    double alpha = mylr.alpha;
    double theta1 = mylr.theta1;
    double theta2 = mylr.theta2;
    double theta3 = mylr.theta3;
    memcpy(out,&alpha,sizeof(alpha));
    out+=sizeof(alpha);
    memcpy(out,&theta1,sizeof(theta1));
    out+=sizeof(theta1);
    memcpy(out,&theta2,sizeof(theta2));
    out+=sizeof(theta2);
    memcpy(out,&theta3,sizeof(theta3));
    out+=sizeof(theta3);
    
    if(tmp_bit>=31){
     out = write_delta_default(in,out,32,length);
    }
    else{
    out = write_delta(delta, out, tmp_bit, length);
    }
    
    
    free(delta);
    
    return out;
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
    //std::cout<<"decompressing all!"<<std::endl;
    double *alpha;
    double *theta1;
    double *theta2;
    double * theta3;
    uint8_t maxerror;
    uint8_t * tmpin=in;
    maxerror = tmpin[0];
    tmpin++;
    alpha = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    theta1 = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    theta2 = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    theta3 = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    if(maxerror>=31){
        read_all_default(tmpin ,0,0, length, maxerror,theta1[0],alpha[0], out);
    }
    else{
        read_all_bit_spline(tmpin ,0,0, length, maxerror,alpha[0],theta1[0],theta2[0],theta3[0], out);
    }

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    double *alpha;
    double *theta1;
    double *theta2;
    double *theta3;
    uint8_t maxerror;
    uint8_t * tmpin=in;
    maxerror = tmpin[0];
    tmpin++;
    alpha = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    theta1 = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    theta2 = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    theta3 = reinterpret_cast<double*>(tmpin);
    tmpin+=8;
    uint32_t tmp=0;
    if(maxerror>=31){
       tmp = read_bit_default(tmpin ,maxerror, l, theta1[0],theta2[0], 0);
     }
    else{
       tmp = read_bit_spline(tmpin ,maxerror, l, alpha[0],theta1[0],theta2[0],theta3[0], 0);
    }

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
uint64_t summation(uint8_t* in, const size_t l, size_t nvalue) {
    return 0;
}
uint32_t get_block_nums(){
      return 1;
}    
std::string name() const {
    return "spine_fix"; 
}    
 void destroy(){} 
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
