
#ifndef PIECEWISERANSAC_H_
#define PIECEWISERANSAC_H_

#include "common.h"
#include "codecs.h"
#include "bit_read.h"
#include "bit_write.h"
#include "lr.h"
#include "RANSAC.h"
#include "rank.h"
#define INF 0x7f7fffff

namespace Codecset {

    
class piecewise_ransac : public IntegerCODEC {
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
  int ransacdelta=20;
void init( int blocks,  int blocksize,int delta){
      block_num=blocks;
      block_size=blocksize;
      temp = ceil((double)block_size/64.);
      ransacdelta = delta;
    
}    
uint8_t* encodeArray8(uint32_t *in, const size_t length,uint8_t*res, size_t nvalue) {
    double *indexes = new double[length];
    double *keys = new double[length];
    uint8_t *out = res;
    uint64_t * writebitmap=new uint64_t[temp];
    for(uint32_t i = 0; i < length; i++){
        indexes[i] = (double) i;
        keys[i] = (double) in[i];
    }
    
    bool *vote = new bool[length];
    //lr.train(indexes, keys, length, 0.0001, 500);
    lr mylr;
    mylr.delta=ransacdelta;
    RANSAC myRansac;
    int usedData = myRansac.compute(mylr,indexes,keys,length,2,vote);
    total_usedData+=usedData;

    //std::cout<<"Theta: "<<mylr.theta0<<" "<<mylr.theta1<<std::endl;
    free(indexes);
    free(keys);
    

    int max_error =0;
    uint64_t tmpbit=0;
    int k=0;
    int wtitebittmp=0;
    int mycounter=0;
    int outlier_num = block_size -usedData ;
    bool have_outlier=true;
    uint32_t *tmpoutlier = new uint32_t[outlier_num];
    if(!outlier_num){
    have_outlier=false;
    }

    
    int *deltax = new int[length];
    int wtiteoutlier=0;
    for(int i=0;i<(long long)length;i++){
        
        if(vote[i]){
            int tmp = (long long) in[i] - (long long)(mylr.theta0+mylr.theta1*(double)mycounter);
            mycounter++;
            deltax[k]=tmp;
            k++;
            if(abs(tmp)>max_error){
                max_error = abs(tmp);
            }
        }
        else{
            tmpbit+= ((1L)<<(63-i%64));
            deltax[k]= 1 ;
            k++;
            
            tmpoutlier[wtiteoutlier]=(long long)in[i];
            wtiteoutlier++;
            
        }
        //if(nvalue==0){
        //std::cout<<"i "<<i<<" vote "<<!vote[i]<<" tmpbit "<<tmpbit<<std::endl;
        //}
        if (i%64 == 63){
        writebitmap[wtitebittmp]= tmpbit;
        wtitebittmp++;
        tmpbit = 0;
        }
    }
    if(block_size%64!=0){
        writebitmap[wtitebittmp]= tmpbit;
        wtitebittmp++;
    }
    int* rank_lut_ = new int [temp];
    int basic_block_size_=64;
    initRankLut(block_size,writebitmap,rank_lut_,basic_block_size_);

    /*
    if(nvalue==0){
    for(int i=0;i<length;i++){
    std::cout<<"i"<<i<<" deltatrue is "<<delta[i]<<std::endl;
    }}
    */
    
    
    /*
        std::map<int,uint32_t>::iterator iter=tmpmap.begin();
        while(iter!=tmpmap.end()){
        std::cout<<iter->first<<" "<<iter->second<<std::endl;
        iter++;
        }
     */
    
    int tmp_bit = bits(max_error)+1;
    
    out[0]=(uint8_t)tmp_bit;
    out++;
    double theta0 = mylr.theta0;
    double theta1 = mylr.theta1;
    memcpy(out,&theta0,sizeof(theta0));
    out+=sizeof(theta0);
    memcpy(out,&theta1,sizeof(theta1));
    out+=sizeof(theta1);
    memcpy(out,&outlier_num,sizeof(outlier_num));
    out+=sizeof(outlier_num);
    memcpy(out,&basic_block_size_,sizeof(basic_block_size_));
    out+=sizeof(basic_block_size_);
    uint8_t *outlier_pos = out;
    out+=sizeof(basic_block_size_);
    
    memcpy(out,writebitmap,temp*sizeof(writebitmap[0]));

    out+=sizeof(writebitmap[0])*temp;
    memcpy(out,rank_lut_,temp*sizeof(rank_lut_[0]));
    out+=sizeof(rank_lut_[0])*temp;
    // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
    

    /*
    if(in[0]==0){
    for(int i=0;i<length;i++){
        std::cout<<i<<" "<<deltax[i]<<" "<<in[i]<<std::endl;
    }
    }
    */
    //std::cout<<"bit_length: "<<tmp_bit<<std::endl;
    
    if(tmp_bit>=31){
        out = write_delta_default(in,out,32,length);
    }
    else{
        out = write_delta(deltax, out, tmp_bit, length);
    }
    
    //delete[] deltax;
    
    int  shift = out - res;
    //std::cout<<" delta length: "<<out - delta_mark<<std::endl;
    
    memcpy(outlier_pos,&shift,sizeof(shift));
    if(have_outlier){
        memcpy(out,tmpoutlier,outlier_num * sizeof(tmpoutlier[0]));

        out+=outlier_num * sizeof(tmpoutlier[0]);
    }
    

    delete[] rank_lut_;
    delete[] deltax;
    delete[] tmpoutlier;
    delete[] vote;
    delete[] writebitmap;
    return out;
    
}

uint32_t *decodeArray8( uint8_t *in, const size_t length, uint32_t *out, size_t nvalue) {
    //std::cout<<"decompressing all!"<<std::endl;
    // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
    double *theta0;
    double *theta1;
    uint8_t maxerror;
    int *outlier_num = 0;
    uint8_t * tmpin=in;
    maxerror = tmpin[0];
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
        /*
        for(int i=0;i<temp;i++){
            std::cout<<bitmap[nvalue][i]<<std::endl;
            
        }
        std::cout<<"outlier"<<std::endl;
        for(int i=0;i<outlier[nvalue].size();i++){
            std::cout<<outlier[nvalue][i]<<std::endl;
        }
        */
        /*
        if(outlier_num==0){
    
        }
        else{uint8_t * outlier_pos = mark + outlier_position;}
        */
        uint8_t * outlier_pos = in + outlier_position[0];
        read_all_bit_ransac(tmpin ,0,0, length, maxerror,theta1[0],theta0[0], out,outlier_pos,outlier_num[0],bitmap_pos);
        //read_all_bit_fix(in ,0,0, length, maxerror[nvalue],mylr.theta1,mylr.theta0, out);
    }

    return out;
}
uint32_t randomdecodeArray8( uint8_t *in, const size_t l, uint32_t *out, size_t nvalue){
    
    //double start =getNow();
    // bit + theta0 + theta1 + #outlier + lookupsize64 + outlier_pos + bitmap + lookup + delta + outlier
    double *theta0;
    double *theta1;
    uint8_t maxerror;
    int *basic_block_size_=0;

    uint8_t * tmpin = in;
    maxerror = tmpin[0];
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
    /*
    uint64_t * bitmap = new uint64_t[temp];
    memcpy(bitmap,bitmap_pos,temp*8);

    int * lookup = new int[temp];
    memcpy(lookup,lookup_pos,temp*4);
    */
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
      //int rankval = rank(l,bitmap,lookup,basic_block_size_);
      //std::cout<<"l "<<l<<" type "<<((tempbitmap>>(63-l%64))&1)<<" rankval "<<rankval<<std::endl;
      if(((tempbitmap[0]>>(63-l%64))&1)){
          uint32_t * tmppoint = reinterpret_cast<uint32_t*>(outlier_pos + (rankval-1) * 4);
          return tmppoint[0];
          
      }
       else{
          tmp = read_bit_ransac(tmpin ,maxerror , l, theta1[0],theta0[0], 0,rankval);
          return tmp;
       }
       
    }
    //double end2 = getNow();
    //std::cout<<"find lower bound time: "<<(end-start)*1000000000<<" read bit time: "<<(end2-start2)*1000000000<<" rate is: "<<(end-start)/(end2-start2)<<std::endl;
    
    
   

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

return total_usedData;
}    
std::string name() const {
    return "piecewise_ransac"; 
}    
void destroy(){}  
};



} // namespace FastPFor

#endif /* SIMDFASTPFOR_H_ */
