


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
#include <rank.h>
#include <popcount.h>

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


long long sum_all_default(uint8_t *in ,int start_byte,int start_index, int numbers,int l){
    long long sum =0;
    for(int i=0;i<numbers;i++){
    uint32_t code = in[0];
    in++;
    code += (in[0]<<8);
    in++;
    code += (in[0]<<16);
    in++;
    code += (in[0]<<24);
    in++;
    sum+=code;
    }
    return sum;
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
    
uint32_t* read_all_bit_double(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out) {

    int left = 0;
    uint64_t decode = 0;
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
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp&((1L<<(l-1))-1);
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l);
        /*
        if(start_index==0 && writeind<=1315){
        std::cout<<"writeind "<<writeind<<" theta0 "<<start_key<<" theta1 "<<slope<<" delta "<<tmpval<<" predict "<<(long long)(start_key + ((double)writeind*slope))<<std::endl;
        }
        */
        tmpval+= (long long)(start_key + ((double)writeind*slope));

        out[writeind]=tmpval;
        writeind++;
        left-=l;
      }
      decode+=((long long)in[start]<<left);
      start++;
      left+=8;
      
    }
    return out+numbers;

}

uint32_t* read_all_bit(uint8_t *in ,int start_byte,int start_index, int numbers,int l,float slope,uint32_t start_key, uint32_t *out) {

    int left = 0;
    uint64_t decode = 0;
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
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp&((1L<<(l-1))-1);
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l);

        tmpval+= ((long long)start_key +(long long) ((float)writeind*(float)slope));

        out[writeind]=tmpval;
        writeind++;
        left-=l;
      }
      decode+=((long long)in[start]<<left);
      start++;
      left+=8;
      
    }
    return out+numbers;

}

void read_all_bit_fix32(uint32_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out) {
    int left = 0;
    uint64_t decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    uint32_t mask0 = (1L<<l)-1;
    uint32_t mask1= (1L<<(l-1))-1;
    uint32_t * res = out;
    end = start +(int)(total_bit/32);
    if(total_bit%32!=0){
      end++;
    }

    while(start<=end){
      /*
      if(writeind>= numbers){
        break;
      }
      */
      while(left>=l){
        uint64_t tmp=decode&mask0;
        long  tmpval = tmp&mask1;
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l); 
        
        tmpval+=  ((long long)start_key +(long long)((double)writeind*slope));
        //out[writeind]=tmpval;
        *res = tmpval;
        res++;
        writeind++;
        left-=l;
        if(left==0){decode=0;}
        //std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
      }
      decode+=((long long)in[start]<<left);
      
      start++;
      left+=32;
      
    }
      
}

void read_all_bit_fix(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out) {
    int left = 0;
    uint64_t decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    uint32_t mask0 = (1L<<l)-1;
    uint32_t mask1= (1L<<(l-1))-1;
    end = start +(int)(total_bit/8);
    uint32_t * res = out;
    if(total_bit%8!=0){
      end++;
    }

    while(start<=end){
      /*
      if(writeind>= numbers){
        break;
      }
      */
      while(left>=l){
        uint64_t tmp=decode&mask0;
        long long tmpval = tmp & mask1;
        if(!(((tmp>>(l-1)) &1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l); 
        tmpval+= (long long)(start_key+(double)writeind*slope);
        *res = tmpval;
        res++;
        writeind++;
        left-=l;
        if(left==0){decode=0;}
        //std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
      }
      decode+=((long long)in[start]<<left);

      
      start++;
      left+=8;
      
    }
      
}


long long sum_all_bit_fix(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope) {
    long long sum=0;
    int left = 0;
    uint64_t decode = 0;
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
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp &((1L<<(l-1))-1);
        if(!(((tmp>>(l-1)) &1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l); 
        sum += tmpval;
        sum += (long long)((double)writeind*slope);
        writeind++;
        left-=l;
        if(left==0){decode=0;}
      }
      decode+=((long long)in[start]<<left);
      
      start++;
      left+=8;
      
    }
    return  sum;
}

void read_all_bit_only(uint8_t *in , int numbers,int l, int *out) {
    int left = 0;
    uint64_t decode = 0;
    int start = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    int end = (int)(total_bit/8);
    if(total_bit%8!=0){
      end++;
    }

    while(start<=end){
      if(writeind>= numbers){
        break;
      }
      while(left>=l){
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp &((1L<<(l-1))-1);
        
        if(!(((tmp>>(l-1)) &1))){
          tmpval=-tmpval;
        }
        decode = (decode>>l); 
        out[writeind]=tmpval;

        writeind++;
        left-=l;
        if(left==0){decode=0;}
        //std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
      }
      decode+=((long long)in[start]<<left);
      
      start++;
      left+=8;
      
    }
      
}

void read_all_bit_nonlinear(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double alpha,double theta1,double theta2, uint32_t *out) {
    int left = 0;
    uint64_t decode = 0;
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
        
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp &((1L<<(l-1))-1);
          
        //std::cout<<"tmp "<<tmp<<" tmpval "<<tmpval<<std::endl;
        if(!(((tmp>>(l-1)) &1))){
          tmpval=-tmpval;
        }
        //std::cout<<"writeind: "<<writeind<<" decode "<<decode<<" tmpval: "<<tmpval<<" l: "<<l<<std::endl;
        decode = (decode>>l); 
       
        tmpval+= (long long) (alpha+theta1*(double)writeind+theta2*(double)writeind*(double)writeind);
        
        out[writeind+start_index]=tmpval;

        writeind++;
        left-=l;
        if(left==0){decode=0;}
        //std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
      }
      decode+=((long long)in[start]<<left);
      
      //std::cout<<"in["<<start<<"] "<<unsigned(in[start])<<std::endl;
      //std::cout<<"decode"<<decode<<std::endl;
      
      start++;
      left+=8;
      
    }
      
}

    
void read_all_bit_spline(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double alpha,double theta1,double theta2,double theta3, uint32_t *out) {
    int left = 0;
    uint64_t decode = 0;
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
        
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp &((1L<<(l-1))-1);
          
        //std::cout<<"tmp "<<tmp<<" tmpval "<<tmpval<<std::endl;
        if(!(((tmp>>(l-1)) &1))){
          tmpval=-tmpval;
        }
        //std::cout<<"writeind: "<<writeind<<" decode "<<decode<<" tmpval: "<<tmpval<<" l: "<<l<<std::endl;
        decode = (decode>>l); 
        tmpval+= (long long) (alpha+theta1*(double)writeind+theta2*(double)writeind*(double)writeind+theta3*(double)writeind*(double)writeind*(double)writeind);
        
        out[writeind+start_index]=tmpval;

        writeind++;
        left-=l;
        if(left==0){decode=0;}
        //std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
      }
      decode+=((long long)in[start]<<left);
      
      //std::cout<<"in["<<start<<"] "<<unsigned(in[start])<<std::endl;
      //std::cout<<"decode"<<decode<<std::endl;
      
      start++;
      left+=8;
      
    }
      
}
    
void read_all_bit_ransac(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out,uint8_t* outlier_pos,int outlier_num, uint8_t* bitmap_pos) {
    int temp = ceil((double)numbers/64.);
    uint64_t * bitmap = new uint64_t[temp];
    memcpy(bitmap,bitmap_pos,temp*8);
    /*
    for(int i=0;i<temp;i++){
    std::cout<<"bitmap "<<i<<" "<<bitmap[i]<<std::endl;
    }
    */
    uint32_t * outlier = new uint32_t[outlier_num];
    memcpy(outlier,outlier_pos,outlier_num*4);
    int left = 0;
    uint64_t decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    end = start +(int)(total_bit/8);
    if(total_bit%8!=0){
      end++;
    }
    int ones=0;
    while(start<=end){
      if(writeind>= numbers){
        break;
      }
      while(left>=l){
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp&((1L<<(l-1))-1);   
        //std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        //std::cout<<"writeind "<<writeind<<" delta is "<<tmpval<<std::endl;
        decode = (decode>>l); 
        //int ranktmp = popcountLinear(bitmap,0,writeind+1);
        
        
        if(((bitmap[writeind/64]>>(63-writeind%64))&1)){
            ones++;
            //std::cout<<"this is a outlier "<<writeind<<"temp rank is "<<ranktmp<<" value is "<<outlier[ranktmp-1]<<std::endl;
            out[writeind+start_index]=outlier[ones-1];
            writeind++;
            
        }
        else{
        tmpval+= (long long) (start_key +((double)(writeind-ones)*slope));
        //std::cout<<"this is not a outlier "<<writeind<<" value is "<<tmpval<<std::endl;
        out[writeind+start_index]=tmpval;
        writeind++;
        }
        
        
        left-=l;
      }
      decode+=((long long)in[start]<<left);
      start++;
      left+=8;
      
    }
   
}       

void read_all_bit_outlier_detection(uint8_t *in ,int start_byte,int start_index, int numbers,int l,double slope,double start_key, uint32_t *out,uint8_t* outlier_pos,int outlier_num, uint8_t* bitmap_pos) {
    int temp = ceil((double)numbers/64.);
    uint64_t * bitmap = new uint64_t[temp];
    memcpy(bitmap,bitmap_pos,temp*8);
    /*
    for(int i=0;i<temp;i++){
    std::cout<<"bitmap "<<i<<" "<<bitmap[i]<<std::endl;
    }
    */
    uint32_t * outlier = new uint32_t[outlier_num];
    memcpy(outlier,outlier_pos,outlier_num*4);
    int left = 0;
    uint64_t decode = 0;
    int start = start_byte;
    int end = 0;
    int total_bit = l*numbers;
    int writeind = 0;
    end = start +(int)(total_bit/8);
    if(total_bit%8!=0){
      end++;
    }
    int ones=0;
    while(start<=end){
      if(writeind>= numbers){
        break;
      }
      while(left>=l){
        uint64_t tmp=decode&((1L<<l)-1);
        long long tmpval = tmp&((1L<<(l-1))-1);   
        //std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
        if(!(tmp>>(l-1))){
          tmpval=-tmpval;
        }
        //std::cout<<"writeind "<<writeind<<" delta is "<<tmpval<<std::endl;
        decode = (decode>>l); 
        //int ranktmp = popcountLinear(bitmap,0,writeind+1);
        
        
        if(((bitmap[writeind/64]>>(63-writeind%64))&1)){
            ones++;
            //std::cout<<"this is a outlier "<<writeind<<"temp rank is "<<ranktmp<<" value is "<<outlier[ranktmp-1]<<std::endl;
            out[writeind+start_index]=outlier[ones-1];
            writeind++;
            
        }
        else{
        tmpval+= (long long) (start_key +((double)(writeind)*slope));
        //std::cout<<"this is not a outlier "<<writeind<<" value is "<<tmpval<<std::endl;
        out[writeind+start_index]=tmpval;
        writeind++;
        }
        
        
        left-=l;
      }
      decode+=((long long)in[start]<<left);
      start++;
      left+=8;
      
    }
   
}       
    
uint32_t read_bit_double(uint8_t *in ,int l ,int to_find, double slope,double start_key,int start) {
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

    uint32_t out = value  +(long long) (start_key+(double)(to_find)*slope);
    return out;      
}
    
uint32_t read_bit(uint8_t *in ,int l ,int to_find, float slope,uint32_t start_key,int start) {
    int start_byte = start+to_find*l/8;
    int start_bit = to_find*l%8;
    int occupy = start_bit;
    uint32_t decode =0;
    int total = 0;
    
    while(total<l){
        uint8_t val = in[start_byte];
        decode +=( (val>>occupy)<<total);
        //std::cout<<"l "<<to_find<<" val "<<val<<" decode "<<decode<<" occupy "<<occupy<<" total "<<total<<std::endl;
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
    int find_bit = to_find*l;
    int start_byte = start+find_bit/8;
    int start_bit = find_bit%8;
    int occupy=start_bit;
    int decode =0;
    int total = 0;
    while(total<l){
        uint8_t val = in[start_byte];
        decode +=((val>>occupy)<<total);
        total += (8-occupy);
        occupy=0;
        start_byte++;
        
     }
    
    decode = decode &((1L<<l)-1);
    bool sign = (decode & (1L<<(l-1)));
    int value = (decode & ((1L<<(l-1))-1)) ;
    if (!sign){value = -value;} 
    //std::cout<<"l: "<<l<<" value: "<<value<<" predict: "<< (long long) (start_key +((float)to_find*slope))<<std::endl;
    uint32_t out = value + (long long)(start_key+(double)to_find*slope);
    return out;   
}

uint32_t read_bit_nonlinear(uint8_t *in ,int l ,int to_find, double alpha,double theta0,double theta1,int start) {
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
    //uint32_t out = value;
    uint32_t out = value + (long long) (alpha +((double)to_find*theta0)+(double)to_find*(double)to_find*theta1);
    return out;   
}

uint32_t read_bit_spline(uint8_t *in ,int l ,int to_find, double alpha,double theta0,double theta1,double theta2,int start) {
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

    uint32_t out = value + (long long) (alpha +((double)to_find*theta0)+(double)to_find*(double)to_find*theta1+(double)to_find*(double)to_find*(double)to_find*theta2);
    return out;   
}
    
uint32_t read_bit_ransac(uint8_t *in ,int l ,int to_find, double slope,double start_key,int start, int rank) {
    int start_byte = start+to_find*l/8;
    int start_bit = to_find*l%8;
    int occupy=start_bit;
    uint64_t decode =0;
    int total = 0;
    
    while(total<l){
        uint8_t val = in[start_byte];
        //std::cout<<"in["<<start_byte<<"] "<<unsigned(in[start_byte])<<std::endl;
        decode +=( (long long)(val>>occupy)<<total);
        total += (8-occupy);
        occupy=0;
        start_byte++;
        
     }
    
    decode = decode &((1L<<l)-1);
    bool sign = (decode & (1L<<(l-1)));
    uint32_t value = (decode & ((1L<<(l-1))-1));
    if (!sign){value = -value;}
    //std::cout<<"l: "<<l<<" value: "<<value<<" predict: "<< (long long) (start_key +((float)to_find*slope))<<std::endl;
    uint32_t out = value + (long long) (start_key +((double)(to_find-rank)*slope));
    return out;   
}

uint32_t read_bit_outlier_detection(uint8_t *in ,int l ,int to_find, double slope,double start_key,int start, int rank) {
    int start_byte = start+to_find*l/8;
    int start_bit = to_find*l%8;
    int occupy=start_bit;
    uint64_t decode =0;
    int total = 0;
    
    while(total<l){
        uint8_t val = in[start_byte];
        //std::cout<<"in["<<start_byte<<"] "<<unsigned(in[start_byte])<<std::endl;
        decode +=( (long long)(val>>occupy)<<total);
        total += (8-occupy);
        occupy=0;
        start_byte++;
        
     }
    
    decode = decode &((1L<<l)-1);
    bool sign = (decode & (1L<<(l-1)));
    uint32_t value = (decode & ((1L<<(l-1))-1));
    if (!sign){value = -value;}
    //std::cout<<"l: "<<l<<" value: "<<value<<" predict: "<< (long long) (start_key +((float)to_find*slope))<<std::endl;
    uint32_t out = value + (long long) (start_key +((double)(to_find)*slope));
    return out;   
}


#if defined(__cplusplus)
}
#endif

#endif 




