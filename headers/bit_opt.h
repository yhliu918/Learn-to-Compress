


#ifndef BIT_H_
#define BIT_H_

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
        int tmp=decode&((1U<<l)-1);
        int tmpval = tmp&((1U<<(l-1))-1);
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l);
        /*
        if(start_index ==0){
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


void read_all_bit_fix(uint8_t *in ,int start_byte,int start_index, int numbers,int l,float slope,float start_key, uint32_t *out) {
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
        int tmp=decode&((1U<<l)-1);
        int tmpval = tmp&((1U<<(l-1))-1);
          
        //std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l);
        
        //std::cout<<"ind: "<<writeind<<" l: "<<l<< " predict: "<<(long long) (start_key +(((float)writeind)*slope))<<" delta: "<<tmpval<<std::endl;
        
        

        tmpval+= (long long) (start_key +((float)writeind*slope));
        
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
    int occupy=start_bit;
    int decode =0;
    uint8_t val1 = in[start_byte];
    /*
    if(start == 143750518){
        std::cout<<"****************"<<std::endl;
        std::cout<<to_find<<","<<l<<","<<start_byte<<","<<start_bit<<","<<unsigned(val1)<<std::endl;
    }
    */
    if (8-occupy>=l){
        decode =((val1>>occupy)&((1U<<l) - 1));
    }
    else{
        start_byte++;
        uint8_t val2 = in[start_byte];
        //std::cout<<start_byte<<","<<unsigned(val2)<<std::endl;
        decode = ((val1>>occupy)+ ((val2 & ((1U<<(l-8+occupy))-1))<<(8-occupy)));

    }
    /*
    if(start == 143750518){
        std::cout<<occupy<<","<<unsigned(decode)<<std::endl;    
    }
    */
    bool sign = (decode & (1U<<(l-1)));
    int value = (decode & ((1U<<(l-1))-1));
    if (!sign){value = -value;}
    uint32_t out = value + start_key +(long long) ((float)(to_find)*slope);
    return out;      
}

#if defined(__cplusplus)
}
#endif

#endif 




