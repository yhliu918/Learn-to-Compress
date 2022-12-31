


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
#include "print128.h"



//given a bit number l(how does it save),should return a vector of numbers
uint8_t *write_delta_s(int *in, uint8_t *out, uint8_t l, int numbers, int size)
{
    uint64_t code = 0;
    int occupy = 0;
    int endbit = (l * numbers);
    int end = 0;
    int *tmpin = in;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = (int)endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    uint64_t left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (tmpin >= in + size)
            {
                occupy = 8;
                break;
            }

            bool sign = 1;
            int tmpnum = tmpin[0];
            if (tmpnum <= 0)
            {
                sign = 0;
                tmpnum = -tmpnum;
            }

            uint64_t value1 = ((tmpnum & ((1L << (l - 1)) - 1)) + (sign << (l - 1)));
            code += (value1 << occupy);
            occupy += l;

            tmpin++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> 8;
            //std::cout<<code<<std::endl;
            code = code & ((1L << 8) - 1);
            occupy -= 8;
            out[0] = unsigned((uint8_t)code);
            code = left_val;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    return out;
}

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



uint8_t * write_delta_T(int *in, uint8_t *out, uint8_t l, int numbers)
{
    uint64_t code = 0;
    int occupy = 0;
    uint64_t endbit = (l * (uint64_t)numbers);
    uint32_t end = 0;
    int writeind = 0;
    int *tmpin = in;
    int readind = 0;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    uint64_t left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (tmpin >= in + numbers)
            {
                occupy = 8;
                break;
            }

            bool sign = 1;
            int tmpnum = tmpin[0];
            if (tmpnum <= 0){
                sign = 0;
                tmpnum=-tmpnum;
            }
            uint64_t value1 =
                (tmpnum & (((1UL) << (uint8_t)(l - 1)) - 1)) 
               + (((uint64_t)sign) << (uint8_t)(l - 1));


            code += (value1 << (uint8_t)occupy);
            occupy += l;
            tmpin++;
            readind++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> (uint8_t)8;
            //std::cout<<code<<std::endl;
            code = code & ((1 << 8) - 1);
            uint8_t tmp_char = code;
            occupy -= 8;
            out[0] = tmp_char;
            code = left_val;
            //std::cout<< writeind<<std::endl;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    
    int pad = ceil((sizeof(uint32_t)*8 - l)/8);
    for (int i = 0; i < pad; i++)
    {
        out[0] = 0;
        out++;
    }
    return out;
}


template <typename T>
uint8_t * write_delta_int_T(std::vector<T>& in,std::vector<bool>& signvec, uint8_t *out, uint8_t l, int numbers)
{
    uint128_t code = 0;
    int occupy = 0;
    uint64_t endbit = (l * (uint64_t)numbers);
    uint64_t end = 0;
    int writeind = 0;
    
    int readind = 0;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    uint64_t left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (readind >= numbers)
            {
                occupy = 8;
                break;
            }

            
            T tmpnum = in[readind];
            bool sign = signvec[readind];
            T value1 =
                (tmpnum & (((T)1 << (uint8_t)(l - 1)) - 1)) 
               + (((T)sign) << (uint8_t)(l - 1));


            code += ((uint128_t)value1 << (uint8_t)occupy);
            occupy += l;
           
            readind++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> (uint8_t)8;
            //std::cout<<code<<std::endl;
            code = code & ((1 << 8) - 1);
            uint8_t tmp_char = code;
            occupy -= 8;
            out[0] = tmp_char;
            code = left_val;
            //std::cout<< writeind<<std::endl;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    
    int pad = 8 - end % 8;
    for (int i = 0; i < pad; i++)
    {
        out[0] = 0;
        out++;
    }
    return out;
}


template <typename T>
uint8_t * write_FOR_int_T(T *in, uint8_t *out, uint8_t l, int numbers)
{
    uint128_t code = 0;
    int occupy = 0;
    uint64_t endbit = (l * (uint64_t)numbers);
    uint64_t end = 0;
    int writeind = 0;
    T *tmpin = in;
    int readind = 0;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    uint64_t left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (tmpin >= in + numbers)
            {
                occupy = 8;
                break;
            }

            
            T tmpnum = tmpin[0];
            T value1 =
                (tmpnum & (((T)1 << l) - 1));


            code += ((uint128_t)value1 << (uint8_t)occupy);
            occupy += l;
            tmpin++;
            readind++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> (uint8_t)8;
            //std::cout<<code<<std::endl;
            code = code & ((1 << 8) - 1);
            uint8_t tmp_char = code;
            occupy -= 8;
            out[0] = tmp_char;
            code = left_val;
            //std::cout<< writeind<<std::endl;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    
    int pad = ceil((sizeof(uint32_t)*8 - l)/8);
    for (int i = 0; i < pad; i++)
    {
        out[0] = 0;
        out++;
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

uint8_t *write_string_delta_string(long_int *in, uint8_t *out, uint32_t l, int numbers)
{
    long_int code = 0;
    int occupy = 0;
    uint64_t endbit = (l * numbers);
    int end = 0;
    long_int *tmpin = in;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = (int)endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    long_int left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (tmpin >= in + numbers)
            {
                occupy = 8;
                break;
            }

            bool sign = 1;
            long_int tmpnum = tmpin[0];
            if (tmpnum <= 0)
            {
                sign = 0;
                tmpnum = -tmpnum;
            }
            long_int value1 = ((tmpnum & (((long_int)1 << (l - 1)) - 1)) + ((long_int)sign << (l - 1)));
            //std::cout<<value1<<std::endl;
            code += (value1 << occupy);
            occupy += l;

            tmpin++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> 8;
            //std::cout<<code<<std::endl;
            code = code & (((long_int)1 << 8) - 1);
            uint8_t tmp_char = code.convert_to<uint8_t>();
            occupy -= 8;
            out[0] = tmp_char;
            code = left_val;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    return out;
}


uint8_t *write_string_delta_string_mpz(std::vector<mpz_t *> in, uint8_t *out, uint32_t l, int numbers)
{
    mpz_t code;
    mpz_init(code);

    int occupy = 0;
    uint64_t endbit = (l * numbers);
    int end = 0;
    int idx = 0;


    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = (int)endbit / 8 + 1;
    }
    uint8_t *last = out + end;


    while (out <= last)
    {
        while (occupy < 8)
        {
            if (idx >=  numbers)
            {
                occupy = 8;
                break;
            }

            bool sign = 1;
            // tmpnum = tmpin[0];

            
            if (mpz_sgn(*in[idx]) <= 0)
            {
                sign = 0;
                mpz_neg(*in[idx], *in[idx]);
            }

            mpz_t value1;
            mpz_init(value1);
            // long_int value1 = ((tmpnum & ((1 << (l - 1)) - 1)) + (sign << (l - 1)));

            mpz_t tmp1;
            mpz_init(tmp1);
            mpz_t tmp2;
            mpz_init(tmp2);
            mpz_set_ui(tmp2, 1);
            mpz_mul_2exp(tmp2, tmp2, l - 1);
            mpz_sub_ui(tmp2, tmp2, 1);
            mpz_and(tmp1, *in[idx], tmp2);
            if (sign)
            {
                mpz_add(value1, tmp1, tmp2);
                mpz_add_ui(value1, value1, 1);

            }
            else
            {
                mpz_set(value1, tmp1);
            }
            mpz_mul_2exp(value1, value1, occupy);
            mpz_add(code, code, value1);

            occupy += l;

            idx++;


        } //end while
        while (occupy >= 8)
        {
            // left_val = code >> 8;
            mpz_t imm_256;
            mpz_init(imm_256);
            mpz_set_ui(imm_256, 256);

            mpz_t tmp_char;
            mpz_init(tmp_char);
            mpz_mod(tmp_char, code, imm_256);
            uint8_t tmp_char_uint8 = mpz_get_ui(tmp_char);
            out[0] = tmp_char_uint8;
            mpz_t tmp_char_left;
            mpz_init(tmp_char_left);
            mpz_div(tmp_char_left, code, imm_256);
            mpz_set(code, tmp_char_left);
            occupy -= 8;
            out++;
            mpz_clear(tmp_char);
            mpz_clear(tmp_char_left);

        }

    }
    return out;
}

uint8_t *write_string_delta_string_128(int128_t *in, uint8_t *out, uint32_t l, int numbers)
{
    uint128_t code = 0;
    int occupy = 0;
    uint64_t endbit = (l * numbers);
    int end = 0;
    int writeind = 0;
    int128_t *tmpin = in;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = (int)endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    uint128_t left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (tmpin >= in + numbers)
            {
                occupy = 8;
                break;
            }

            bool sign = 1;
            int128_t tmpnum = tmpin[0];
            if (tmpnum <= 0)
            {
                sign = 0;
                tmpnum = -tmpnum;
            }
            uint128_t value1 = ((tmpnum & (( ((uint128_t)1) << (l - 1)) - 1)) + ( ((uint128_t)sign) << (l - 1)));
           
            code += (value1 << occupy);
            occupy += l;

            tmpin++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> 8;
            //std::cout<<code<<std::endl;
            code = code & ((1 << 8) - 1);
            uint8_t tmp_char = code;
            occupy -= 8;
            out[0] = tmp_char;
            code = left_val;
            writeind ++;
            //std::cout<< writeind<<std::endl;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    return out;
}

template <typename T>
uint8_t * write_delta_string(T *in, std::vector<bool>& signvec, uint8_t * string_len, uint8_t *out, uint32_t l, size_t numbers)
{
    T code = 0;
    int occupy = 0;
    uint64_t endbit = ((l+8) * numbers);
    int end = 0;
    int writeind = 0;
    T *tmpin = in;
    int readind = 0;
    if (endbit % 8 == 0)
    {
        end = endbit / 8;
    }
    else
    {
        end = (int)endbit / 8 + 1;
    }
    uint8_t *last = out + end;
    T left_val = 0;

    while (out <= last)
    {
        while (occupy < 8)
        {
            if (tmpin >= in + numbers)
            {
                occupy = 8;
                break;
            }

            bool sign = signvec[readind];
            T tmpnum = tmpin[0];
            T value1 =
                (tmpnum & ((((T)1) << (uint8_t)(l - 1)) - 1)) 
               + (((T)sign) << (uint8_t)(l - 1));
            value1 += ((T)string_len[readind])<<(uint8_t)l;


            code += (value1 << (uint8_t)occupy);
            occupy += (l+8);
            tmpin++;
            readind++;
        } //end while
        while (occupy >= 8)
        {
            left_val = code >> (uint8_t)8;
            //std::cout<<code<<std::endl;
            code = code & ((1 << 8) - 1);
            uint8_t tmp_char = code;
            occupy -= 8;
            out[0] = tmp_char;
            code = left_val;
            //std::cout<< writeind<<std::endl;
            //std::cout<<occupy<<" "<<left_val<<" "<<unsigned(out[0])<<std::endl;
            out++;
        }
    }
    
    int pad = ceil((sizeof(T)*8 - l)/8);
    for (int i = 0; i < pad; i++)
    {
        out[0] = 0;
        out++;
    }
    return out;
}






#endif
