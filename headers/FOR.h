/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */
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


  
uint32_t * encodeArray(uint32_t *in, const size_t length,uint32_t *out, size_t &nvalue) {
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
                                      uint32_t *out, size_t &nvalue) {
    nvalue = in[0];
    ++in;
    if(nvalue == 0) return in;
    uint32_t m = in[0];
    ++in;
    uint32_t M = in[0];
    ++in;
    int b = bits(static_cast<uint32_t>(M-m));
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(uint32_t k = 0; k<nvalue/32; ++k) {
        unpack32[b](m,in+b*k,out+32*k);
    }
    out = out + nvalue/32*32;
    in = in + nvalue/32*b;

    for(uint32_t k=nvalue/32*32; k+16<=nvalue; k+=16,out+=16) {
        in = unpack16[b](m,in,out);
    }
    for(uint32_t k=nvalue/16*16; k+8<=nvalue; k+=8,out+=8) {
        in = unpack8[b](m,in,out);
    }
    // we could pack the rest, but we don't  bother
    for(uint32_t k=nvalue/8*8; k<nvalue; ++k,in++,out++) {
        out[0] = in [0];
    }
    return in;
}
uint32_t randomdecodeArray( uint32_t *in, const size_t l,
                                      uint32_t *out, size_t &nvalue){
    nvalue = in[0];
    ++in;
    if(nvalue == 0) return in[0];
    uint32_t m = in[0];
    ++in;
    uint32_t M = in[0];
    ++in;
    int b = bits(static_cast<uint32_t>(M-m));
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    //std::cout<<" min: "<<m<<" max: "<<M<<"bit: "<<b<<std::endl;

    uint32_t recover=0;
    in = in + ((int)l/32)*b;
    uint32_t number_left = l - ((int)l/32)*32;
    uint32_t number_occupy = (number_left*b)/32;
    in+= number_occupy;
    if(b==32){
        return in[0]+m;
    }
    
    uint32_t bit_left = number_left*b - number_occupy*32;
    //std::cout<<"number_left: "<<number_left<<" number_occupy: "<<number_occupy<<" bit_left: "<<bit_left<<std::endl;
    
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

}
    
uint8_t* encodeArray8( uint32_t *in, const size_t length, uint8_t *out,
                   size_t &nvalue) {
    std::cout<<"Haven't implement. Please try uint32_t one..."<<std::endl;
    return out;
}
uint32_t *decodeArray8( uint8_t *in, const size_t length,
                              uint32_t *out, size_t &nvalue) {
    std::cout<<"Haven't implement. Please try uint32_t one..."<<std::endl;
    return out;
}
uint32_t randomdecodeArray8(uint8_t *in, const size_t l,uint32_t *out, size_t &nvalue){
    std::cout<<"Haven't implement. Please try uint32_t one..."<<std::endl;
    return 1;
}  
std::string name() const {
    return "FrameofReference"; 
}    
  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
