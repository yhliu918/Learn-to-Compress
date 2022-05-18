
#ifndef ELIASFANO_H_
#define ELIASFANO_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#define INF 0x7f7fffff

namespace Codecset {
    
class eliasfana : public IntegerCODEC {
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
  //double * models;
  
void init(int blocks, int blocksize,int extra){
      block_num=blocks;
      block_size=blocksize;
      //models = new double [block_num*2];
}

int random(int m){
        return rand()%m;
}

// bit + theta0 + theta1 + delta   
uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {
    int compressedBufferPointer2 = 0;

    uint8_t * out = res;
    uint64_t buffer1 = 0;
    int bufferLength1 = 0;
    uint64_t buffer2 = 0;
    int bufferLength2 = 0;
    uint32_t largestblockID = in[length - 1];
    double averageDelta = (double)largestblockID / (double)length;

    double averageDeltaLog = log(averageDelta);
    int lowBitsLength = averageDeltaLog < 0 ? 0 : averageDeltaLog;
    uint64_t lowBitsMask = (((uint64_t)1 << lowBitsLength) - 1);
    int compressedBufferPointer1 = 0;
    compressedBufferPointer2 = lowBitsLength * length/8 + 6;

    memcpy(out, &length, sizeof(int));
    out += sizeof(int);
    memcpy(out, &lowBitsLength, sizeof(uint8_t));
    out += sizeof(uint8_t);
    int lastDocID = 0;
    for(int i=0;i<length;i++){
        // docID strictly monotone/increasing numbers, docIdDelta strictly positive, no zero allowed
        int docID = in[i];
        uint docIdDelta = (docID - lastDocID - 1); 

        // low bits
        // Store the lower l= log(u / n) bits explicitly
        // binary packing/bit packing

        buffer1 <<= lowBitsLength;
        buffer1 |= (docIdDelta & lowBitsMask);
        bufferLength1 += lowBitsLength;

        // flush buffer to compressedBuffer
        while (bufferLength1 > 7)
        {
            bufferLength1 -= 8;
            memcpy(out, &(buffer1 >> bufferLength1), sizeof(uint8_t));    
            out += sizeof(uint8_t);    
        }

 

        // high bits
        // Store high bits as a sequence of unary coded gaps
        // 0=1, 1=01, 2=001, 3=0001, ...
        // https://en.wikipedia.org/wiki/Unary_coding

        // length of unary code 
        uint32_t unaryCodeLength = (uint32_t)(docIdDelta >> lowBitsLength) + 1; 
        buffer2 <<= (int)unaryCodeLength;
        
        // set most right bit 
        buffer2 |= 1;
        bufferLength2 += (int)unaryCodeLength;

        // flush buffer to compressedBuffer
        while (bufferLength2 > 7)
        {
            bufferLength2 -= 8; 
            memcpy(out, &(buffer2 >> bufferLength2), sizeof(uint8_t));
            out += sizeof(uint8_t);
        }

        lastDocID = docID;
    }
    if (bufferLength1 > 0)
    {
        memcpy(out, &(buffer1 << (8 - bufferLength1)), sizeof(uint8_t));
        out += sizeof(uint8_t);
    }

    if (bufferLength2 > 0)
    {
        memcpy(out, &(buffer2 << (8 - bufferLength2)), sizeof(uint8_t));
        out += sizeof(uint8_t);
    }

    std::cout<<"delta: "<<averageDelta<<"\n";
    std::cout<<"lowBitsLength: "<<lowBitsLength<<"\n";
    std::cout<<"bits/DocID"<<(double)compressedBufferPointer2 * (double)8 / (double)length <<" ("<< (2+averageDeltaLog)<<")\n";
    std::cout<<"compressed rate "<<(double)compressedBufferPointer2/((double)length* 4)<<"\n";
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
    //std::cout<<"decompressing all!"<<std::endl;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    

}
/*
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
*/

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
    return "eliasfano"; 
}    
 void destroy(){} 
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
