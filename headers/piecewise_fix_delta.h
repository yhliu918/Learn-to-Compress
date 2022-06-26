
#ifndef PIECEWISEFIX_DELTA_H_
#define PIECEWISEFIX_DELTA_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {




    
class piecewise_fix_delta : public IntegerCODEC {
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

// bit + theta0 + theta1 + delta[0]+ delta(after delta)
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {
    uint8_t *out = res;

    std::vector<double> indexes;
    std::vector<double> keys;
    for(uint32_t i = 0; i < length; i++){
        indexes.emplace_back((double) i);
        keys.emplace_back((double) in[i]);
    }
    int *delta = new int[length];
    
    lr mylr;
    mylr.caltheta(indexes,keys,length);

    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - (long long)(mylr.theta0+mylr.theta1*(double)i);
        delta[i]=tmp;
        
    }
    
    int *delta2 = new int[length];
    uint32_t max_error =abs(delta[0]); //delta[0], delta[1]-delta[0],....
    delta2[0]=delta[0];
    for(int i=1;i<(long long) length;i++){
        delta2[i]=delta[i]-delta[i-1]; 
        if(abs(delta2[i])>max_error){
            max_error = abs(delta2[i]);
        }
    }
    int tmp_bit = bits(max_error)+1;
    if(max_error==0){
            tmp_bit=0;
    }
    //std::cout<<"tmp_bit you need is "<<tmp_bit<<std::endl;

    out[0]=(uint8_t)tmp_bit;
    out++;
    double theta0 = mylr.theta0;
    double theta1 = mylr.theta1;
    memcpy(out,&theta0,sizeof(theta0));
    out+=sizeof(theta0);
    memcpy(out,&theta1,sizeof(theta1));
    out+=sizeof(theta1);
    free(delta);

    if(tmp_bit==0){
        free(delta2);
        return out;
    }

    else if(tmp_bit>=31){
        out = write_delta_default(in,out,32,length);
    }
    else{
        out = write_delta(delta2, out, tmp_bit, length);
    }
    

    free(delta2);
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
    if(maxerror==0){
        for(uint32_t i=0;i<length;i++){
            out[i] = (long long) (theta0 +((double)i*theta1));
        }
        return out;
    }
    if(maxerror>=31){
        read_all_default(tmpin ,0,0, length, maxerror,theta1,theta0, out);
        return out;
    }
    else{
        int * delta = (int*)malloc(block_size * sizeof(int));
        read_all_bit_only(tmpin ,length, maxerror, delta);
        for(int i=1;i<block_size;i++){
            delta[i] = delta[i-1]+delta[i];
        }
        for(int i=0;i<block_size;i++){
            out[i]=delta[i]+(long long) (theta0 +((double)i*theta1));
        }
        delete [] delta;
        return out;
    }
    
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    
    std::cout<<"Piecewise_fix_delta doesn't support random access..."<<std::endl;
    return 0;

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
    return "piecewise_fix_delta"; 
}    
 void destroy(){} 
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
