


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
void read_all_bit(uint8_t *in ,int start_byte,int start_index, int numbers,int l, int *delta) {
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
        delta[writeind+start_index]=tmpval;

        writeind++;
        left-=l;
      }
      decode+=(in[start]<<left);
      start++;
      left+=8;
      
    }
    

}


int read_bit(uint8_t *in ,int l ,int to_find, int start) {
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
        /*
        if(start == 143750518){
            std::cout<<start_byte<<","<<unsigned(val2)<<std::endl;  
        }
        */
    }
    /*
    if(start == 143750518){
        std::cout<<occupy<<","<<unsigned(decode)<<std::endl;    
    }
    */
    bool sign = (decode & (1U<<(l-1)));
    int value = (decode & ((1U<<(l-1))-1));
    if (!sign){value = -value;}
    return value;      
}

#if defined(__cplusplus)
};
#endif

#endif 




