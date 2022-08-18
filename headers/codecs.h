
#ifndef CODECS_H_
#define CODECS_H_

#include "common.h"
#include "util.h"


namespace Codecset {

class NotEnoughStorage : public std::runtime_error {
public:
  size_t required; // number of 32-bit symbols required
  NotEnoughStorage(const size_t req)
      : runtime_error(""), required(req) {}
};

class IntegerCODEC {
public:
  /**
   * You specify input and input length, as well as
   * output and output length. nvalue gets modified to
   * reflect how much was used. If the new value of
   * nvalue is more than the original value, we can
   * consider this a buffer overrun.
   *
   * You are responsible for allocating the memory (length
   * for *in and nvalue for *out).
   */
  std::vector<int> mul_add_diff_set;
  std::vector<int> mul_add_diff_set_minus;
  virtual void init(int blocks, int blocksize, int extra)=0;
  virtual uint32_t * encodeArray( uint32_t *in, const size_t length,
                           uint32_t *out, size_t nvalue) = 0;
  virtual uint8_t * encodeArray8( uint32_t *in, const size_t length,
                           uint8_t *out,  size_t nvalue) = 0;
  /**
   * Usage is similar to decodeArray except that it returns a pointer
   * incremented from in. In theory it should be in+length. If the
   * returned pointer is less than in+length, then this generally means
   * that the decompression is not finished (some scheme compress
   * the bulk of the data one way, and they then they compress remaining
   * integers using another scheme).
   *
   * As with encodeArray, you need to have length element allocated
   * for *in and at least nvalue elements allocated for out. The value
   * of the variable nvalue gets updated with the number actually use
   * (if nvalue exceeds the original value, there might be a buffer
   * overrun).
   */
  virtual  uint32_t *decodeArray( uint32_t *in, const size_t length,
                                      uint32_t *out, size_t nvalue) = 0;
  virtual  uint32_t randomdecodeArray( uint32_t *in, const size_t l,
                                      uint32_t *out, size_t nvalue) = 0;    
  virtual  uint32_t *decodeArray8( uint8_t *in, const size_t length,
                                      uint32_t *out, size_t nvalue) = 0;
  virtual  uint32_t randomdecodeArray8( uint8_t *in, const size_t l,
                                      uint32_t *out, size_t nvalue) = 0;  
  virtual  uint64_t summation( uint8_t *in, const size_t l, size_t nvalue) = 0;  
  virtual  uint32_t get_block_nums() = 0;  
  virtual  void destroy() = 0;  
  virtual ~IntegerCODEC() {}
    


  virtual std::string name() const = 0;
};

/******************
 * This just copies the data, no compression.
 */
class JustCopy : public IntegerCODEC {
public:
    

  int block_num;
  int block_size;
  
void init(int blocks, int blocksize,int extra){
      block_num=blocks;
      block_size=blocksize;
}
  uint32_t* encodeArray( uint32_t *in, const size_t length, uint32_t *out,
                   size_t nvalue) {
    memcpy(out, in, sizeof(uint32_t) * length);
    nvalue = length;
    return out;
  }
  uint8_t* encodeArray8( uint32_t *in, const size_t length, uint8_t *out,
                   size_t nvalue) {
    std::cout<<"Haven't implement. Please try uint32_t one..."<<std::endl;
    return out;
  }
  // like encodeArray, but we don't actually copy
  void fakeencodeArray( uint32_t * /*in*/, const size_t length,
                       size_t nvalue) {
    nvalue = length;
  }

   uint32_t *decodeArray( uint32_t *in, const size_t length,
                              uint32_t *out, size_t nvalue) {
    memcpy(out, in, sizeof(uint32_t) * length);
    nvalue = length;
    return in + length;
  }
  uint32_t *decodeArray8( uint8_t *in, const size_t length,
                              uint32_t *out, size_t nvalue) {
    std::cout<<"Haven't implement. Please try uint32_t one..."<<std::endl;
    return out;
  }
  uint32_t randomdecodeArray(uint32_t *in, const size_t l,uint32_t *out, size_t nvalue){
      std::cout<<"not implement yet!"<<std::endl;
      return 1;
  }
  uint32_t randomdecodeArray8(uint8_t *in, const size_t l,uint32_t *out, size_t nvalue){
    std::cout<<"Haven't implement. Please try uint32_t one..."<<std::endl;
    return 1;
  }
  
  uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    std::cout<<"Haven't implement yet..."<<std::endl;
    return 0;
  }
  
  uint32_t get_block_nums(){
      return 1;
  }
  void destroy(){}
  std::string name() const { return "JustCopy"; }
};



} // namespace FastPFor

#endif /* CODECS_H_ */
