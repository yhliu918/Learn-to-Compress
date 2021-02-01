


#ifndef BIT_READ_H_
#define BIT_READ_H_

#include <stdio.h>
#include<iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <cmath>
#include <time.h>
#include<algorithm>
#include<queue>
#include<vector>
#include<cmath>
#include <getopt.h>


#ifdef __cplusplus
extern "C" {
#endif
//given a bit number l(how does it save),should return a vector of numbers
    

void read_all_default(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out){
    for(int i=0;i<numbers;i++){
    uint32_t code = in[0];
    in++;
    code += (in[0]<<8);
    in++;
    code += (in[0]<<16);
    in++;
    code += (in[0]<<24);
    in++;
    out[0]=code;
    out++;
    }
}
    
uint32_t read_bit_default(uint8_t *in ,int l ,int to_find, double slope,double start_key,int start) {
    in+=4*to_find;
    uint32_t code = in[0];
    in++;
    code += (in[0]<<8);
    in++;
    code += (in[0]<<16);
    in++;
    code += (in[0]<<24);
    in++;
    
    return code;   
}
    
void read_all_bit(uint8_t *in ,int start_byte,int start_index, int numbers,int l,float slope,uint32_t start_key, uint32_t *out) {
    int left = 0;
    int decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    end = start +(int)(total_bit/8);
    if(total_bit%8!=0){
      end++;
    }

    while(start<=end){
      if(writeind>= numbers){
        break;
      }
      while(left>=l){
        int tmp=decode&((1L<<l)-1);
        long long tmpval = tmp&((1L<<(l-1))-1);
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l);
        /*
        if(start_key==3875996905){
        std::cout<<"slope: "<<slope<<" start key: "<<start_key<< " predict: "<<(long long)start_key +(long long) ((float)writeind * slope)<<" delta: "<<tmpval<<std::endl;
       }
        */
        tmpval+= (start_key +(long long) ((float)writeind*slope));
        
        out[writeind+start_index]=tmpval;

        writeind++;
        left-=l;
      }
      decode+=(in[start]<<left);
      start++;
      left+=8;
      
    }
    

}


void read_all_bit_fix(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out) {
    int left = 0;
    int decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    end = start +(int)(total_bit/8);
    if(total_bit%8!=0){
      end++;
    }

    while(start<=end){
      if(writeind>= numbers){
        break;
      }
      while(left>=l){
        int tmp=decode&((1L<<l)-1);
        int tmpval = tmp&((1L<<(l-1))-1);
          
        //std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l); 

        tmpval+= (long long) (start_key +((double)writeind*slope));
        
        out[writeind+start_index]=tmpval;

        writeind++;
        left-=l;
      }
      decode+=(in[start]<<left);
      start++;
      left+=8;
      
    }
      
}

    
void read_all_bit_ransac(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out,std::map<int,uint32_t> outlier) {
    int left = 0;
    int decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    end = start +(int)(total_bit/8);
    if(total_bit%8!=0){
      end++;
    }

    while(start<=end){
      if(writeind>= numbers){
        break;
      }
      while(left>=l){
        int tmp=decode&((1L<<l)-1);
        int tmpval = tmp&((1L<<(l-1))-1);
          
        //std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l); 
        std::map<int,uint32_t>::iterator iter;
        
        iter = outlier.find(writeind);
        while(iter!=outlier.end()){
            out[writeind+start_index]=iter->second;
            writeind++;
            iter = outlier.find(writeind);
        }

        tmpval+= (long long) (start_key +((double)writeind*slope));
        
        out[writeind+start_index]=tmpval;

        writeind++;
        left-=l;
      }
      decode+=(in[start]<<left);
      start++;
      left+=8;
      
    }
    
   
}
    
uint32_t read_bit(uint8_t *in ,int l ,int to_find, float slope,uint32_t start_key,int start) {
    int start_byte = start+to_find*l/8;
    int start_bit = to_find*l%8;
    int occupy = start_bit;
    int decode =0;
    int total = 0;
    
    while(total<l){
        uint8_t val = in[start_byte];
        decode +=( (val>>occupy)<<total);
        total += (8-occupy);
        occupy=0;
        start_byte++;
        
     }
    
    decode = decode &((1L<<l)-1);
    bool sign = (decode & (1L<<(l-1)));
    int value = (decode & ((1L<<(l-1))-1));
    if (!sign){value = -value;}
    //std::cout<<"l: "<<l<<" value: "<<value<<" predict: "<< start_key +(long long) ((float)(to_find)*slope)<<std::endl;
    uint32_t out = value + start_key +(long long) ((float)(to_find)*slope);
    return out;      
}

uint32_t read_bit_fix(uint8_t *in ,int l ,int to_find, double slope,double start_key,int start) {
    int start_byte = start+to_find*l/8;
    int start_bit = to_find*l%8;
    int occupy=start_bit;
    int decode =0;
    int total = 0;
    
    while(total<l){
        uint8_t val = in[start_byte];
        decode +=( (val>>occupy)<<total);
        total += (8-occupy);
        occupy=0;
        start_byte++;
        
     }
    
    decode = decode &((1L<<l)-1);
    bool sign = (decode & (1L<<(l-1)));
    int value = (decode & ((1L<<(l-1))-1));
    if (!sign){value = -value;}
    //std::cout<<"l: "<<l<<" value: "<<value<<" predict: "<< (long long) (start_key +((float)to_find*slope))<<std::endl;
    uint32_t out = value + (long long) (start_key +((double)to_find*slope));
    return out;   
}
#if defined(__cplusplus)
}
#endif

#endif 




