
#ifndef PIECEWISEFIX32_H_
#define PIECEWISEFIX32_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {




    
class piecewise_fix_32 : public IntegerCODEC {
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

int random(int m){
        return rand()%m;
}

// bit + theta0 + theta1 + delta   
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];
    //double *keys_sample = new double [length];
    //double *indexes_sample = new double[length];
    uint8_t *out = res;
    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    int *delta = new int[length];
    /*
    keys_smooth[0]=keys[0];
    for(uint32_t i = 1; i < length-1; i++){
        double tmpsum =0;
        for(int j=-1;j<=1;j++){
            if(i+j<0){
                tmpsum+=keys[0];
            }
            else if(i+j>=length){
                tmpsum+=keys[length-1];
            }
            else{
                tmpsum+=keys[i+j];
            }
        }
        keys_smooth[i]=(double)tmpsum/3.;
    }
    keys_smooth[length-1]=keys[length-1];
    */
    /*
    int sample_length=0;
    for(int i=0;i<length;i++){
        if(random(3)==0){
            indexes_sample[sample_length]=i;
            keys_sample[sample_length]=keys[i];
            sample_length++;
        }

    }
    */
    //lr.train(indexes, keys, length, 0.0001, 500);

    lr mylr;
    //mylr.caltheta_LOO(indexes,keys,length);
    mylr.caltheta(indexes,keys,length);
    
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
    double theta0 = mylr.theta0;
    double theta1 = mylr.theta1;
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
       tmp = read_bit_default(tmpin ,maxerror, l, theta1,theta0, 0);
     }
    else{
       tmp = read_bit_fix(tmpin ,maxerror, l, theta1,theta0, 0);
    }

    return tmp;
    
   

}
uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    long long sum =0;
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
        sum = sum_all_default(tmpin ,0,0, nvalue, maxerror);
    }
    else{
        //uint32_t total = nvalue *(nvalue-1)/2;
        sum = sum_all_bit_fix(tmpin ,0,0, nvalue, maxerror,theta1);
        sum +=((long long) theta0 )*nvalue;
        //sum +=(long long)(theta1*total);
    }
    return (uint64_t)sum;
  }

uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *res,
                   size_t nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];
    uint32_t *out = res;
    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    int *delta = new int[length];

    lr mylr;
    mylr.caltheta(indexes,keys,length);
    
    free(indexes);
    free(keys);
    int max_error =0;
    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - ((long long)mylr.theta0+(long long)(mylr.theta1*(double)i));
        delta[i]=tmp;
        if(abs(tmp)>max_error){
            max_error = abs(tmp);
        }
    }

    int tmp_bit = bits(max_error)+1;
    out[0]=tmp_bit;
    out++;
    double theta0 = mylr.theta0;
    double theta1 = mylr.theta1;
    memcpy(out,&theta0,sizeof(theta0));
    out+=2;
    memcpy(out,&theta1,sizeof(theta1));
    out+=2;
    
    if(tmp_bit>=31){
     for(int i=0;i<block_size;i++){
            out[0]=in[i];
            out++;
        }
        free(delta);
        return out;
    }
    else{
    out = write_delta32(delta, out, tmp_bit, length);
    }
    
    
    free(delta);
    
    return out;
}
uint32_t *decodeArray( uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue) {
    double theta0;
    double theta1;
    int maxerror;
    uint32_t * tmpin=in;
    maxerror = tmpin[0];
    tmpin++;
    memcpy(&theta0,tmpin,8);
    tmpin+=2;
    memcpy(&theta1,tmpin,8);
    tmpin+=2;
    if(maxerror>=31){
        memcpy(out,in,sizeof(uint32_t)*length);
        return out;
    }
    else{
        read_all_bit_fix32(tmpin ,0,0, length, maxerror,theta1,theta0, out);
    }

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
    return "piecewise_fix_32"; 
}    
 void destroy(){} 
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
