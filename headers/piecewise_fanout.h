
#ifndef PIECEWISE_FANOUT_H_
#define PIECEWISE_FANOUT_H_

#include "common.h"
#include "codecs.h"
#include "time.h"
#include "bit_read.h"
#include "bit_write.h"
#include "caltime.h"
#include "split.h"
#define INF 0x7f7fffff

namespace Codecset {

bool greater( const seg &v1, const seg &v2)
{
    return v1.start < v2.start;
}


    
class piecewise_fanout : public IntegerCODEC {
public:
  using IntegerCODEC::encodeArray;
  using IntegerCODEC::decodeArray;
  using IntegerCODEC::randomdecodeArray;
  using IntegerCODEC::encodeArray8;
  using IntegerCODEC::decodeArray8;
  using IntegerCODEC::randomdecodeArray8;
  using IntegerCODEC::init;
   
  std::vector<uint8_t*> block_start_vec;
  std::vector<uint32_t> segment_index;
  uint32_t total_byte =0;
  int block_num;
  int block_size;

//start_index + bit + theta0 + theta1 + numbers + delta
void init( int blocks,  int blocksize,int delta){
      block_num=blocks;
      block_size=blocksize;
    
}   


uint32_t lower_bound( uint32_t v,uint32_t len)
{
    uint32_t m;
    uint32_t x=0;
    uint32_t y=len-1;
    while(x <= y)
    {
    
        m = x+(y-x)/2;
        if(v<segment_index[m]) y = m-1;
        else x = m+1;
    }
    return y;
   
}

uint8_t * encodeArray8(uint32_t *in, const size_t length,uint8_t *res, size_t nvalue) {

  double *indexes = new double[length];
  double *keys = new double[length];
  for(int i = 0; i < length; i++){
      indexes[i] = (double) i  ;
      keys[i] = (double) in[i];
  }
  std::queue<seg> q;
  std::vector<seg> v;
    
  seg initseg;
  initseg.start = 0;
  initseg.end = length-1;
  initseg.caltheta(indexes,keys,length);
  initseg.calbit(indexes,keys,length);
  q.push(initseg);
    
  while(!q.empty()){
      seg tmpseg = q.front();
      
      int mid = ( tmpseg.start + tmpseg.end )/2;
      seg subseg1;
      subseg1.start = tmpseg.start;
      subseg1.end = mid;
      subseg1.caltheta(indexes+(subseg1.start),keys+(subseg1.start),(subseg1.end-subseg1.start+1));
      subseg1.calbit(indexes+(subseg1.start),keys+(subseg1.start),(subseg1.end-subseg1.start+1));
      
      seg subseg2;
      subseg2.start = mid+1;
      subseg2.end = tmpseg.end;
      subseg2.caltheta(indexes+(subseg2.start),keys+(subseg2.start),(subseg2.end-subseg2.start+1));
      subseg2.calbit(indexes+(subseg2.start),keys+(subseg2.start),(subseg2.end-subseg2.start+1));
      q.pop();
      if(tmpseg.totalbyte > subseg1.totalbyte+subseg2.totalbyte){
          q.push(subseg1);
          q.push(subseg2);
      }
      else{
          v.push_back(tmpseg);
      }
      
  }
  free(indexes);
  free(keys);
  std::sort(v.begin(),v.end(),greater);
  for(int i=0;i<(int)v.size();i++){
      v[i].start += nvalue * block_size;
      v[i].end += nvalue * block_size;
      int numbers = v[i].end - v[i].start+1;
      //std::cout<<cal_score_tmp(in+v[i].start-nvalue * block_size, numbers)<<std::endl;
      
      uint8_t * descriptor = (uint8_t*)malloc(numbers * sizeof(uint64_t)+1024);
      uint8_t *out = descriptor;
            
      memcpy(out,&(v[i].start),sizeof(v[i].start));
      out+=sizeof(v[i].start);
            
      out[0]=(uint8_t)v[i].tmpbit;
      out++;
      
      memcpy(out,&v[i].theta0,sizeof(v[i].theta0));
      out+=sizeof(v[i].theta0);
      memcpy(out,&v[i].theta1,sizeof(v[i].theta1));
      out+=sizeof(v[i].theta1);
      memcpy(out,&numbers,sizeof(numbers));
      out+=sizeof(numbers);

      if(v[i].tmpbit>=32){
         
         out = write_delta_default(in+v[i].start-nvalue*length,out,32,numbers);
      }
      else{
        out=write_delta(v[i].delta, out, v[i].tmpbit, numbers);
      }
      
      free(v[i].delta);
      
      total_byte += (out-descriptor);
      segment_index.push_back(v[i].start);
      descriptor = (uint8_t*)realloc(descriptor,(out-descriptor));
      block_start_vec.push_back(descriptor);
  }
  
    
    return res;
    
}
    
uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
//start_index + bit + theta0 + theta1 + numbers + delta   
    uint32_t * tmpout = out;
    for(int i=0;i<(int)block_start_vec.size();i++){
        uint8_t * this_block = block_start_vec[i];
        uint8_t * tmpin=this_block;
        
        int start_ind;
        double theta0;
        double theta1;
        uint8_t maxerror;
        int numbers;
        
        memcpy(&start_ind,tmpin,4);
        tmpin+=4;
        memcpy(&maxerror,tmpin,1);
        tmpin++;
        memcpy(&theta0,tmpin,sizeof(theta0));
        tmpin+=sizeof(theta0);
        memcpy(&theta1,tmpin,sizeof(theta1));
        tmpin+=sizeof(theta1);
        memcpy(&numbers,tmpin,4);
        tmpin +=4;
        
            if(maxerror>=32){
                
                read_all_default(tmpin ,0,0, numbers, maxerror,theta1,theta0, tmpout);
                tmpout+=numbers;
            }
            else{
                tmpout = read_all_bit_double(tmpin ,0,start_ind,numbers,maxerror,theta1,theta0, tmpout);
            }

    }
    
        
      
    

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){

    uint32_t length=segment_index.size();
    uint8_t * this_block = block_start_vec[lower_bound(l,length)];

    uint8_t * tmpin = this_block;
    int start_ind;
    double theta0;
    double theta1;
    uint8_t maxerror;
    int numbers;
     
    memcpy(&start_ind,tmpin,4);
    tmpin+=4;
    memcpy(&maxerror,tmpin,1);
    tmpin++;
    memcpy(&theta0,tmpin,sizeof(theta0));
    tmpin+=sizeof(theta0);
    memcpy(&theta1,tmpin,sizeof(theta1));
    tmpin+=sizeof(theta1);
    memcpy(&numbers,tmpin,4);
    tmpin +=4;

    //std::cout<<theta0<<" , "<<theta1<<std::endl;

    //std::cout<< "indexing "<<l<<std::endl;
    uint32_t tmp;
    if(maxerror>=32){
         tmp = read_bit_default(tmpin ,maxerror , l-start_ind,theta1,theta0,0);
    }
    else{
         
         tmp = read_bit_double(tmpin ,maxerror , l-start_ind,theta1,theta0,0);
    }
    return tmp;

}
uint64_t summation( uint8_t *in, const size_t l, size_t nvalue){
    
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
   
    std::cout<<"Total block num is "<<block_start_vec.size()<<std::endl;
    
     return total_byte;
 }    
void destroy(){
for(int i=0;i<(int)block_start_vec.size();i++){
       free(block_start_vec[i]);
   }
    
}
std::string name() const {
    return "piecewise_fanout"; 
}    
  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
