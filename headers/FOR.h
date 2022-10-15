
#ifndef FOR_H_
#define FOR_H_

#include "common.h"
#include "codecs.h"
#include "bpacking.h"
#include "forutil.h"


namespace Codecset {


class FOR : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
  using IntegerCODEC::init;
    using IntegerCODEC::summation;

  int block_num;
  int block_size;
  
void init(int blocks, int blocksize,int extra){
      block_num=blocks;
      block_size=blocksize;
}
  
uint32_t * encodeArray(uint32_t *in, const size_t length,uint32_t *res, size_t nvalue) {
    uint32_t * out = res;
    out[0] = length;
    ++out;

    if(length == 0) return out;
    uint32_t m = in[0];
    uint32_t M = in[0];
    for(uint32_t i = 1; i < length; ++i) {
        if(in[i]>M) M=in[i];
        if(in[i]<m) m=in[i];
    }
    int b = bits(static_cast<uint32_t>(M-m));
    if(b==31){
        b=32;
    }
    
    out[0] = m;
    ++out;
    out[0] = M;
    ++out;
    uint32_t k = 0;
    //std::cout<<m<<","<<M<<std::endl;
    for(; k+32<=length; k+=32,in+=32) {
        
        out = pack32[b](m,in,out);
    }
    for(; k+16<=length; k+=16,in+=16) {
        //std::cout<<m<<std::endl;
        out = pack16[b](m,in,out);
    }
    for(; k+8<=length; k+=8,in+=8) {

        out = pack8[b](m,in,out);
    }
    // we could pack the rest, but we don't  bother
    for(; k<length; ++k,in++,out++) {
        out[0] = in [0];
    }
    return out;
    
}
uint32_t *decodeArray( uint32_t *in, const size_t length,
                                      uint32_t *out,  size_t nvalue) {
    uint32_t*res = out;
    nvalue = in[0];
    ++in;
    if(nvalue == 0) return in;
    uint32_t m = in[0];
    ++in;
    uint32_t M = in[0];
    ++in;
    
    int b = bits(static_cast<uint32_t>(M-m));
    if(b==31){
        b=32;
    }
    //std::cout<<"bit "<<b<<" min "<<m<<" max "<<M<<std::endl;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    if(b ==0){
        for(int i=0;i<(int)length;i++){
            out[i]=m;
        }
        return out;
    }
    for(uint32_t k = 0; k<length/32; ++k) {
        unpack32[b](m,in+b*k,res+32*k);
        
    }
    res = res + length/32*32;
    in = in + length/32*b;

    for(uint32_t k=length/32*32; k+16<=length; k+=16,res+=16) {
        in = unpack16[b](m,in,res);
    }
    for(uint32_t k=length/16*16; k+8<=length; k+=8,res+=8) {
        in = unpack8[b](m,in,res);

    }
    // we could pack the rest, but we don't  bother
    for(uint32_t k=length/8*8; k<length; ++k,in++,res++) {
        res[0] = in [0];
    }
    return res;
}
uint32_t randomdecodeArray( uint32_t *in, const size_t l,
                                      uint32_t *out,  size_t nvalue){
    uint32_t tmpval = in[0];
    ++in;
    if(tmpval == 0) return in[0];
    uint32_t m = in[0];
    ++in;
    uint32_t M = in[0];
    ++in;
    int b = bits(static_cast<uint32_t>(M-m));
    if(b==0){
        return m;
    }
    if(b==31){
        b=32;
    }
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    
    uint32_t recover=0;
    in = in + ((int)l/32)*b;
    int number_left = l - ((int)l/32)*32;
    //std::cout<<"nvalue "<<nvalue<<" number_left "<<number_left<<" nvalue - l "<<nvalue - l<<std::endl;
    //std::cout<<"to_find "<<l<<" block_size "<<block_size<<" bit "<<b<<std::endl;
    
    
    // if(((int)l/32) == ((int)nvalue/32)){

    //     /*
    //     if( number_left>=16  ){
    //         number_left -= 16;
    //         in = in +(int)ceil((double)b*16./32.);
    //         if(number_left>=8){
    //             number_left -= 8;
    //             in = in +(int)ceil((double)b*8./32.);
    //         }
    //         else{
    //             unpack8[b](m,in,out);
    //             return out[number_left];
    //         }
        
    //     }
    //     else if(number_left>=8){
    //         number_left -= 8;
    //         in = in +(int)ceil((double)b*8./32.);
    //         return in[number_left];
            
    
    //     }
    
    //     if(number_left>0){
    //         return in[number_left];
    //     }
    //     */
    //     uint32_t *res = new uint32_t[32];
    //     uint32_t *tmpres =res;
    //     for(uint32_t k=nvalue/32*32; k+16<=nvalue; k+=16,res+=16) {
    //         in = unpack16[b](m,in,res);
    //     }
    //     for(uint32_t k=nvalue/16*16; k+8<=nvalue; k+=8,res+=8) {
    //         in = unpack8[b](m,in,res);
    //     }
    // // we could pack the rest, but we don't  bother
    //     for(uint32_t k=nvalue/8*8; k<nvalue; ++k,in++,res++) {
    //         res[0] = in [0];
    //     }
    //     recover = tmpres[number_left];
    //     free(tmpres);
    //     return recover;
    // }
    // else{
        uint32_t number_occupy = (number_left*b)/32;
        in+= number_occupy;

        if(b==32){ 
            return in[0]+m;
        }
    

        long long bit_left = number_left*b - number_occupy*32;

        if(32-bit_left>=b){
            recover = (in[0]>>bit_left) & ((1U<<b)-1) ;
            recover += m ;
        
            return recover;
        }
        else{
            recover = (( in[1]&(((1U<<(b+bit_left-32))-1)) )<<(32-bit_left)) + (in[0]>>bit_left);
            recover +=m;
            return recover;
        }
    
    
    // }

 /*
    uint32_t number_occupy = (number_left*b)/32;
    in+= number_occupy;
    if(b==32){ 
        return in[0]+m;
    }
    

    long long bit_left = number_left*b - number_occupy*32;
    
    //std::cout<<"ind "<<l<<" min: "<<m<<" max: "<<M<<"bit: "<<b<<std::endl;
    //std::cout<<"number_left: "<<number_left<<" number_occupy: "<<number_occupy<<" bit_left: "<<bit_left<<std::endl;
    //std::cout<<"to the end of this seg "<<block_size - l-1<<std::endl;

    if(32-bit_left>=b){
        recover = (in[0]>>bit_left) & ((1U<<b)-1) ;
        recover += m ;
        
        return recover;
    }
    else{
        recover = (( in[1]&(((1U<<(b+bit_left-32))-1)) )<<(32-bit_left)) + (in[0]>>bit_left);
        recover +=m;
        return recover;
    }
*/
}
    
uint8_t* encodeArray8( uint32_t *in, const size_t length, uint8_t *res,
                    size_t nvalue) {
    uint32_t * out=reinterpret_cast<uint32_t*>(res);
    uint32_t * mark_out = out;
    out = encodeArray(in,length,out,nvalue);
    
    res = reinterpret_cast<uint8_t*>(mark_out);
    uint8_t *tmp_res = reinterpret_cast<uint8_t*>(out);
    return tmp_res;
}
uint32_t *decodeArray8( uint8_t *in, const size_t length,
                              uint32_t *out,  size_t nvalue) {
    uint32_t * tmpin = reinterpret_cast<uint32_t*>(in);
    return decodeArray(tmpin,length,out,nvalue);

}
uint32_t randomdecodeArray8(uint8_t *in, const size_t l,uint32_t *out, size_t nvalue){
    uint32_t * tmpin = reinterpret_cast<uint32_t*>(in);
    uint32_t tmp = randomdecodeArray(tmpin,l,out,nvalue);
    //std::cout<<tmp<<std::endl;
    return tmp;
} 
uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    uint32_t* res = (uint32_t*)malloc(nvalue*sizeof(uint32_t));
    decodeArray8(in,nvalue,res,nvalue);
    uint64_t sum=0;
    for(int i=0;i<(int)nvalue;i++){
        sum+=res[i];
    }
    return sum;
}
uint32_t get_block_nums(){
      return 1;
}
void destroy(){}
std::string name() const {
    return "FrameofReference"; 
}    
  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
