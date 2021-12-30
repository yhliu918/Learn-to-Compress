
#ifndef RANSACFIX_H_
#define RANSACFIX_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#include "RANSAC.h"
#include "rank.h"
#define INF 0x7f7fffff

namespace Codecset {

    
class ransac_fix : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
  using IntegerCODEC::init;

  int temp;
  int total_usedData=0;
  int block_num;
  int block_size;
  int ransacdelta=20;
void init( int blocks,  int blocksize,int delta){
      block_num=blocks;
      block_size=blocksize;
      temp = ceil((double)block_size/64.);
      ransacdelta = delta;
    
}    
uint8_t* encodeArray8(uint32_t *in, const size_t length,uint8_t*res, size_t nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];
    uint8_t *out = res;
    uint64_t * writebitmap=new uint64_t[temp];
    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    int *delta = new int[length];
    bool *vote = new bool[length];
    //lr.train(indexes, keys, length, 0.0001, 500);
    lr mylr;
    mylr.delta=ransacdelta;
    RANSAC myRansac;
    int usedData = myRansac.compute(mylr,indexes,keys,length,2,vote);
    total_usedData+=usedData;

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

    int tmp_bit = bits(max_error)+1;
    out[0]=(uint8_t)tmp_bit;
    out++;
    double theta0 = mylr.theta0;
    double theta1 = mylr.theta1;
    //models[nvalue*2] = theta0;
    //models[nvalue*2+1] = theta1;
    memcpy(out,&theta0,sizeof(theta0));
    out+=sizeof(theta0);
    memcpy(out,&theta1,sizeof(theta1));
    out+=sizeof(theta1);

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
    double theta0;
    double theta1;
    uint8_t maxerror;
    uint8_t * tmpin=in;
    memcpy(&maxerror,tmpin,1);
    tmpin++;
    memcpy(&theta0,tmpin,8);
    tmpin+=8;
    memcpy(&theta1,tmpin,8);
    tmpin+=8;
    if(maxerror>=31){
        read_all_default(tmpin ,0,0, length, maxerror,theta1,theta0, out);
    }
    else{
        
        read_all_bit_fix(tmpin ,0,0, length, maxerror,theta1,theta0, out);
    }

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    double theta0;
    double theta1;
    uint8_t maxerror;
    uint8_t * tmpin=in;
    memcpy(&maxerror,tmpin,1);
    tmpin++;
    memcpy(&theta0,tmpin,8);
    tmpin+=8;
    memcpy(&theta1,tmpin,8);
    tmpin+=8;
    uint32_t tmp=0;
    if(maxerror>=31){
        //uint32_t * interpret = reinterpret_cast<uint32_t*>(tmpin);
       tmp = read_bit_default(tmpin ,maxerror, l, theta1,theta0, 0);
       //return interpret[l];
     }
    else{
       tmp = read_bit_fix(tmpin ,maxerror, l, theta1,theta0, 0);
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
uint32_t get_block_nums(){

return total_usedData;
}    
std::string name() const {
    return "ransac_fix"; 
}    
void destroy(){}  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
