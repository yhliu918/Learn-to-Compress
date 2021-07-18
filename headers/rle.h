
#ifndef RLE_H_
#define RLE_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {


  
class rle : public IntegerCODEC {
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



// key + times
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {
    uint32_t * out=reinterpret_cast<uint32_t*>(res);
    uint32_t * mark_out = out;
    out = encodeArray(in,length,out,nvalue);
    res = reinterpret_cast<uint8_t*>(mark_out);
    uint8_t *tmp_res = reinterpret_cast<uint8_t*>(out);
    return tmp_res;
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
    uint32_t * tmpin = reinterpret_cast<uint32_t*>(in);
    return decodeArray(tmpin,length,out,nvalue);
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    
    uint32_t * tmpin = reinterpret_cast<uint32_t*>(in);

    uint32_t tmp = 0;
    tmp = randomdecodeArray(tmpin,l,out,nvalue);
    
    return tmp;

}

uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *res,
                   size_t nvalue) {

    uint32_t *out = res+1;
    uint32_t total_pair=0;
    uint32_t lastkey = in[0];
    uint32_t count = 0;
    out[0]=lastkey;
    out++;

    for(uint32_t i =0;i<length;i++){
        if(in[i]==lastkey){
            count++;
        }
        else{
            out[0]=count;
            out++;
            total_pair ++;
            
            
            lastkey = in[i];
            out[0]=lastkey;
            out++;
            count = 1;
        }
    }
   out[0]=count;
   out++;
   total_pair ++;

   res [0] = total_pair;


   return out;
}
uint32_t *decodeArray( uint32_t *in, const size_t length,
                              uint32_t *res, size_t nvalue) {
    uint32_t * tmpin =in;
    uint32_t * out =res;
    uint32_t total_pair=tmpin[0];
    tmpin++;
    for(uint32_t i=0;i<total_pair;i++){
        uint32_t key = tmpin[0];
        
        tmpin++;
        uint32_t count = tmpin[0];
        tmpin++;
        
        for(uint32_t j =0;j<count;j++){
            out[0]=key;
            out++;
        }
    }
    /*
    if(nvalue ==617){
        for(int i=0;i<length;i++){
        std::cout<<"out["<<i<<"] is "<<res[i]<<std::endl;
        }
    }
    */
    return res;
}
uint32_t randomdecodeArray(uint32_t *in, const size_t l,uint32_t *out, size_t nvalue){
    uint32_t * tmpin =in;
    uint32_t tmp = 0;
    uint32_t total_pair=tmpin[0];
    tmpin++;
    uint32_t totalcount =0;
    for(uint32_t i=0;i<total_pair;i++){
        uint32_t key = tmpin[0];
        tmpin++;
        uint32_t count = tmpin[0];
        tmpin++;
        totalcount+=count;
        if(totalcount>l){
            tmp = key;
            break;
        }
        
    }    
    return tmp;
}
uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    return 0;
}
uint32_t get_block_nums(){
      return 1;
}    
std::string name() const {
    return "rle"; 
}    
 void destroy(){} 
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
