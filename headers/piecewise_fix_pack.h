
#ifndef PIECEWISEFIX_PACK_H_
#define PIECEWISEFIX_PACK_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {




    
class piecewise_fix_pack : public IntegerCODEC {
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

    uint8_t *out = res;
    std::vector<double> indexes;
    std::vector<double> keys;
    for(uint32_t i = 0; i < length; i++){
        indexes.emplace_back((double) i);
        keys.emplace_back((double) in[i]);
    }
    int *delta = new int[length];
    lr mylr;
    //mylr.caltheta_LOO(indexes,keys,length);
    mylr.caltheta(indexes,keys,length);
    
    //std::cout<<"Theta: "<<mylr.theta0<<" "<<mylr.theta1<<std::endl;

    int max_error =0;
    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - ((long long)mylr.theta0+(long long)(mylr.theta1*(double)i));
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

uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *res,size_t nvalue) {

    uint32_t *out = res;
    std::vector<double> indexes;
    std::vector<double> keys;
    for(uint32_t i = 0; i < length; i++){
        indexes.emplace_back((double) i);
        keys.emplace_back((double) in[i]);
    }
    int *delta = new int[length];
    lr mylr;
    mylr.caltheta(indexes,keys,length);

    int max_error =0;
    //double theta1_i = 0.0;
    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - (long long)(mylr.theta0+mylr.theta1*(double)i);
        //theta1_i+=mylr.theta1;
        delta[i]=tmp;
        if(abs(tmp)>max_error){
            max_error = abs(tmp);
        }
    }
    int tmp_bit = bits(max_error)+1;
    uint32_t* unsigned_delta = new uint32_t[length];
    uint32_t* mark = unsigned_delta;
    for(int i=0;i<(long long)length;i++){
        if(delta[i]>0){
            unsigned_delta[i]= ((delta[i] & ((1L<<(tmp_bit-1))-1)) + (1<<(tmp_bit-1)));
        }
        else{
            unsigned_delta[i] = -delta[i];
        }
    }
    
    out[0]=tmp_bit;
    out++;
    double theta0 = mylr.theta0;
    double theta1 = mylr.theta1;
    memcpy(out,&theta0,sizeof(theta0));
    out+=2;
    memcpy(out,&theta1,sizeof(theta1));
    out+=2;
    if (tmp_bit>=31){
        for(int i=0;i<block_size;i++){
            out[0]=in[i];
            out++;
        }
        free(delta);
        free(unsigned_delta);
        return out;
    }

    uint32_t k = 0;
    for(; k+32<=length; k+=32,unsigned_delta+=32) {
        out = pack32[tmp_bit](0,unsigned_delta,out);
    }
    for(; k+16<=length; k+=16,unsigned_delta+=16) {
        //std::cout<<m<<std::endl;
        out = pack16[tmp_bit](0,unsigned_delta,out);
    }
    for(; k+8<=length; k+=8,unsigned_delta+=8) {

        out = pack8[tmp_bit](0,unsigned_delta,out);
    }
    // we could pack the rest, but we don't  bother
    for(; k<length; ++k,unsigned_delta++,out++) {
        out[0] = unsigned_delta[0];
    }
    free(delta);
    free(mark);
    return out;
}
uint32_t *decodeArray( uint32_t *in, const size_t length,uint32_t *out, size_t nvalue) {
    double theta0;
    double theta1;
    int maxerror;
    uint32_t * tmpin=in;
    
    memcpy(&maxerror,tmpin,4);
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
        //int *delta = new int[length];
        uint32_t *delta = new uint32_t[length];
        uint32_t * mark = delta;
        for(uint32_t k = 0; k<length/32; ++k) {
            unpack32[maxerror](0,tmpin+maxerror*k,delta+32*k);
        }
        delta = delta + length/32*32;
        tmpin = tmpin + length/32*maxerror;

        for(uint32_t k=length/32*32; k+16<=length; k+=16,delta+=16) {
            tmpin = unpack16[maxerror](0,tmpin,delta);
        }
        for(uint32_t k=length/16*16; k+8<=length; k+=8,delta+=8) {
            tmpin = unpack8[maxerror](0,tmpin,delta);
        }
    // we could pack the rest, but we don't  bother
        for(uint32_t k=length/8*8; k<length; ++k,in++,delta++) {
            delta[0] = tmpin [0];
        }
        
        //double theta1_i =0.0;
        for(uint32_t i=0;i<length;i++){
            out[i] = (long long)(theta0+theta1*(double)i);
            //theta1_i +=theta1;
            /*
            if(!(((mark[i]>>(maxerror-1)) &1))){
                out[i]-=mark[i];
            }
            else{
                out[i]+=(mark[i]-(1<<(maxerror-1)));
            }
            */
            out[i]+=mark[i];
            
        }
        free(mark);
        return out;
    }
   
}
uint32_t randomdecodeArray(uint32_t *in, const size_t l,uint32_t *out, size_t nvalue){
    std::cout<<"Haven't implement. Please try uint8_t one..."<<std::endl;
    return 1;
}
uint32_t get_block_nums(){
      return 1;
}    
std::string name() const {
    return "piecewise_fix_pack"; 
}    
 void destroy(){} 
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
