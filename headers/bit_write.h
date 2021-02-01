


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
    int code =0;
    int occupy = 0;
    int endbit = (l*numbers);
    int end=0;

    if(endbit%8==0){
        end=endbit/8;
    }
    else{
        end = (int)endbit/8+1;
    }
    uint8_t* last=out+end;
    uint32_t left_val = 0;

    while(out<=last){

        
        while(occupy<8){
            
            bool sign = 1;
            if (in[0] <= 0){
                sign = 0;
                in[0]=-in[0];
            }

            uint32_t value1= ((in[0] & ((1L<<(l-1))-1)) + (sign<<(l-1)));
            code += (value1<<occupy);
            occupy += l;

            in++;
            
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




