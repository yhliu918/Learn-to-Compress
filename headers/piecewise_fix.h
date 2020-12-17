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
#include "time.h"
#include "bit_opt.h"
#include "LinearRegression.h"
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
  LinearRegression lr;
  uint8_t maxerror = 0;



uint8_t* write_delta(int *in,uint8_t* out, uint8_t l, int numbers){
    int code =0;
    int occupy = 0;

    uint8_t left_val = 0;

    for(int i = 0; i < numbers;i++){

        code = left_val;
        while(occupy<8){
            
            bool sign = 1;
            if (in[i] <= 0){
                sign = 0;
                in[i]=-in[i];
            }

            uint8_t value1= ((in[i] & ((1U<<(l-1))-1)) + (sign<<(l-1)));
            code += (value1<<occupy);
            occupy += l;

            i++;
            
        }//end while
        if(occupy>8){
            left_val = code >> 8;
            code = code & ((1U<<8) - 1);
        }
        else{
            left_val=0;
        
        }
        
        occupy = occupy % 8;
        
        out[0]= unsigned((uint8_t)code);
        out++;
        if(i==numbers){
            if (left_val){
                out[0]=(uint8_t)left_val;
                out++;
            }
            if(left_val==0&&occupy!=0){
                out[0]=(uint8_t)left_val;
                out++;
            }
            
        }
        i--;
    }

    
    return out;

}
    
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *out, size_t &nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];

    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    int *delta = new int[length];

    lr.train(indexes, keys, length, 0.01, 1500);
    uint8_t max_error =0;
    for(int i=0;i<length;i++){
        int tmp = (int) in[i] - (int)lr.predict(indexes[i]);
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
          
    out = write_delta(delta, out, tmp_bit, length);
    maxerror = max_error;
    
    
    return out;
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t &nvalue) {
    std::cout<<"decompressing all!"<<std::endl;

    int start_byte=0;
    int end_index =0;
    int start_index = 0;
    int index=0;
    read_all_bit(in ,0,0, length, maxerror,lr.theta[1],lr.theta[0], out);

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t &nvalue){

    uint32_t tmp = read_bit(in ,maxerror , l, lr.theta[1],lr.theta[0], 0);
    return tmp;

}

uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *out,
                   size_t &nvalue) {
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return out;
}
uint32_t *decodeArray( uint32_t *in, const size_t length,
                              uint32_t *out, size_t &nvalue) {
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return out;
}
uint32_t randomdecodeArray(uint32_t *in, const size_t l,uint32_t *out, size_t &nvalue){
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return 1;
}
    
std::string name() const {
    return "FrameofReference"; 
}    
  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
