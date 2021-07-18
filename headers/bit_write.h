


#ifndef BIT_WRITE_H_
#define BIT_WRITE_H_

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
    
uint8_t* write_delta(int *in,uint8_t* out, uint8_t l, int numbers){
    uint64_t code =0;
    int occupy = 0;
    int endbit = (l*numbers);
    int end=0;
    int* tmpin =in;
    if(endbit%8==0){
        end=endbit/8;
    }
    else{
        end = (int)endbit/8+1;
    }
    uint8_t* last=out+end;
    uint64_t left_val = 0;

    while(out<=last){
        
        while(occupy<8){
            
            bool sign = 1;
            int tmpnum = tmpin[0];
            if (tmpnum <= 0){
                sign = 0;
                tmpnum=-tmpnum;
            }

            uint64_t value1= ((tmpnum & ((1L<<(l-1))-1)) + (sign<<(l-1)));
            code += (value1<<occupy);
            occupy += l;

            tmpin++;
            
        }//end while
        while(occupy>=8){
            left_val = code >> 8;
            //std::cout<<code<<std::endl;
            code = code & ((1L<<8) - 1);
            occupy-=8;
            out[0]= unsigned((uint8_t)code);
            code = left_val;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    

    
    return out;

}

uint32_t* write_delta32(int *in,uint32_t* out, uint8_t l, int numbers){
    uint64_t code =0;
    int occupy = 0;
    int endbit = (l*numbers);
    int end=0;
    int* tmpin =in;
    if(endbit%32==0){
        end=endbit/32;
    }
    else{
        end = (int)endbit/32+1;
    }
    uint32_t* last=out+end;
    uint64_t left_val = 0;

    while(out<=last){
        
        while(occupy<32){
            
            bool sign = 1;
            int tmpnum = tmpin[0];
            if (tmpnum <= 0){
                sign = 0;
                tmpnum=-tmpnum;
            }

            uint64_t value1= ((tmpnum & ((1L<<(l-1))-1)) + (sign<<(l-1)));
            code += (value1<<occupy);
            occupy += l;

            tmpin++;
            
        }//end while
        while(occupy>=32){
            left_val = code >> 32;
            //std::cout<<code<<std::endl;
            code = code & ((1L<<32) - 1);
            occupy-=32;
            out[0]= (uint32_t)code;
            code = left_val;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    

    
    return out;

}


uint8_t* write_delta_default(uint32_t *in,uint8_t* out, uint8_t l, int numbers){
    for(int i=0;i<numbers;i++){
        uint32_t code =in[i];
        out[0]=code & ((1L<<8)-1);
        out++;
        out[0]=(code>>8) &((1L<<8)-1);
        out++;
        out[0]=(code>>16) &((1L<<8)-1);
        out++;
        out[0]=(code>>24) &((1L<<8)-1);
        out++;
    }

    
    return out;

}

#if defined(__cplusplus)
}
#endif

#endif 




