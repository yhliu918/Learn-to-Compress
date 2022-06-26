
#ifndef PIECEWISE_VARILENGTH_H_
#define PIECEWISE_VARILENGTH_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#include "RANSAC.h"
#include "rank.h"
#define INF 0x7f7fffff

namespace Codecset {

    
class piecewise_varilength : public IntegerCODEC {
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
void init( int blocks,  int blocksize,int delta){
      block_num=blocks;
      block_size=blocksize;
      temp = ceil((double)block_size/64.);
      

    
}    
int cal_log(int x){
    int tmp_bit = bits(x)+1;
    tmp_bit = std::min(32,tmp_bit);
    return tmp_bit;

}
uint8_t* encodeArray8(uint32_t *in, const size_t length,uint8_t*res, size_t nvalue) {
    std::vector<double> indexes;
    std::vector<double> keys;
    for(uint32_t i = 0; i < length; i++){
        indexes.emplace_back((double) i);
        keys.emplace_back((double) in[i]);
    }
    
    
    lr mylr;
    mylr.caltheta(indexes,keys,length);
    
    //std::cout<<"Theta: "<<mylr.theta0<<" "<<mylr.theta1<<std::endl;
    
    int *counter = new int[32];
    int *delta = new int[length];
    int max_error =0;
    for(int i=0;i<32;i++){
        counter[i]=0;
    }
    for(int i=0;i<(long long)length;i++){
        int tmp = (long long) in[i] - (long long)(mylr.theta0+mylr.theta1*(double)i);
        delta[i]=tmp;
        counter[cal_log(abs(tmp))]++;
        if(abs(tmp)>max_error){
            max_error = abs(tmp);
        }
    }
    int max_bit = cal_log(max_error);
    
    int vari_length =0;
    for(int i=0;i<32;i++){
        vari_length += counter[i] * i;
    }
    vari_length+=  (cal_log(max_bit)-1)* block_size;

    int quantile_sum=0;
    //int threshold =0;
    int compress_len =0;
    int compress_min= block_size*32;
    max_bit =0;
    for(int i=0;i<32;i++){
        quantile_sum+=counter[i];
        compress_len = quantile_sum * i + (block_size-quantile_sum)*32+temp*12*8+12*8;
        if(quantile_sum==block_size){
            compress_len = quantile_sum * i;
            max_bit = i;
        }
        if(compress_len<compress_min){
            compress_min = compress_len;
            //threshold=i;
        }
        if(quantile_sum==block_size){
            break ;
        }
    }
    if (vari_length<compress_len){
        total_usedData += (vari_length/8);
        total_usedData+=16;

    }
    else{
        total_usedData +=(compress_len/8) ;
    }
    
    
    

    delete [] counter;
    delete [] delta;

    return res;
    
}

uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
    //std::cout<<"decompressing all!"<<std::endl;
    // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
    
    uint8_t whether_outlier;
    uint8_t * tmpin=in;
    whether_outlier=(tmpin[0]>>7); //1000 0000
    if(whether_outlier){
        double *theta0;
        double *theta1;
        uint8_t maxerror;
        int *outlier_num = 0;
        maxerror = tmpin[0]-(1L<<7);
        tmpin++;
        theta0 = reinterpret_cast<double*>(tmpin);
        tmpin+=8;
        theta1 = reinterpret_cast<double*>(tmpin);
        tmpin+=8;
        outlier_num = reinterpret_cast<int*>(tmpin);
        tmpin +=8;
        int* outlier_position;
        outlier_position = reinterpret_cast<int*>(tmpin);
        tmpin+=4;
        uint8_t * bitmap_pos = tmpin;
        tmpin+=temp*12;
        //std::cout<<"maxerror "<<unsigned(maxerror)<<" theta0 "<<theta0<<" theta1 "<<theta1<<" outlier_num "<<outlier_num<<" outlier_position "<<outlier_position<<std::endl;
    

        if(maxerror>=31){
            read_all_default(tmpin ,0,0, length, maxerror,theta1[0],theta0[0], out);
        }
        else{
            uint8_t * outlier_pos = in + outlier_position[0];
            read_all_bit_outlier_detection(tmpin ,0,0, length, maxerror,theta1[0],theta0[0], out,outlier_pos,outlier_num[0],bitmap_pos);
        //read_all_bit_fix(in ,0,0, length, maxerror[nvalue],mylr.theta1,mylr.theta0, out);
        }

        return out;
    }
    else{
        double theta0;
        double theta1;
        uint8_t maxerror;
        
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
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    
    //double start =getNow();
    // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
    uint8_t * tmpin = in;

    uint8_t whether_outlier;
    whether_outlier=(tmpin[0]>>7); //1000 0000
    if(whether_outlier){
        double *theta0;
        double *theta1;
        uint8_t maxerror;
        int *basic_block_size_=0;
        maxerror = tmpin[0]-(1L<<7);
        tmpin++;
        theta0 = reinterpret_cast<double*>(tmpin);
        tmpin+=8;
        theta1 = reinterpret_cast<double*>(tmpin);
        tmpin+=12;
        basic_block_size_ = reinterpret_cast<int*>(tmpin);
        tmpin+=4;
        int* outlier_position;
        outlier_position = reinterpret_cast<int*>(tmpin);
        tmpin+=4;
        uint8_t * bitmap_pos = tmpin;
        tmpin+=temp*8;
        uint8_t * lookup_pos = tmpin;
        tmpin+=temp*4;
        uint8_t * outlier_pos = in + outlier_position[0];
    
        //std::cout<<"maxerror "<<unsigned(maxerror)<<" theta0 "<<theta0<<" theta1 "<<theta1<<" outlier_num "<<outlier_num<<" outlier_position "<<outlier_position<<std::endl;
        uint32_t tmp=0;
        if(maxerror>=31){
            tmp = read_bit_default(tmpin ,maxerror , l, theta1[0],theta0[0], 0);
            return tmp;
        }
        else{
            //std::cout<<"num "<<l<<" type "<<((bitmap[nvalue][l/64]>>(63-l%64))&1)<<std::endl;
            int fetch_pos = (l/basic_block_size_[0]);
            int * lookup_num ;
            lookup_num = reinterpret_cast<int*>(lookup_pos + 4 * fetch_pos);
            uint64_t * tempbitmap;
            tempbitmap = reinterpret_cast<uint64_t*>(bitmap_pos + 8 * fetch_pos);
            int rankval =lookup_num[0] + popcountLinear(tempbitmap, 0 , (l&(basic_block_size_[0] - 1))+1);
      
            if(((tempbitmap[0]>>(63-l%64))&1)){
                uint32_t * tmppoint = reinterpret_cast<uint32_t*>(outlier_pos + (rankval-1) * 4);
                return tmppoint[0];
          
            }
            else{
                tmp = read_bit_outlier_detection(tmpin ,maxerror , l, theta1[0],theta0[0], 0,rankval);
                return tmp;
            }
       
        }
    

    }
    else{
        double theta0;
        double theta1;
        uint8_t maxerror;
        
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
uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    return 0;
}
uint32_t get_block_nums(){

return total_usedData;
}    
std::string name() const {
    return "piecewise_varilength"; 
}    
void destroy(){}  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
