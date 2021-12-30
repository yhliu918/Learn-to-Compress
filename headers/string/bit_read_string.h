
#include <algorithm>
#include <cmath>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>
#include "string_utils.h"

void read_all_bit_fix_string(uint8_t *in, int start_byte, int start_index, int numbers, int l, long_int slope, long_int start_key, long_int *out)
{

    int left = 0;
    long_int decode = 0;
    uint64_t start = start_byte;
    uint64_t end = 0;
    uint64_t total_bit = l * numbers;
    int writeind = 0;
    end = start + (int)(total_bit / 8);
    long_int *res = out;
    if (total_bit % 8 != 0)
    {
        end++;
    }
    while (start <= end)
    {

        while (left >= l)
        {
            long_int tmp = decode & (((long_int)1 << l) - 1);
            long_int tmpval = tmp & (((long_int)1 << (l - 1)) - 1);
            
            if (!(((tmp >> (l - 1)) & (long_int)1)))
            {
                tmpval = -tmpval;
            }
            decode = (decode >> l);
            tmpval += (start_key + (writeind * slope));
            
            *res = tmpval;
            //std::cout<<"index "<<writeind<<" tmpval "<<res[0]<<std::endl;
            res++;
            
            writeind++;
            left -= l;
            if (left == 0)
            {
                decode = 0;
            }
            //std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
        }
        decode += ((long_int)in[start] << left);

        start++;
        left += 8;
    }
    
}

void read_outlier_detect_string(uint8_t *in, int start_byte, int start_index, int numbers, int l, long_int slope, long_int start_key, long_int *out, uint8_t *outlier_pos, int outlier_num, uint8_t *bitmap_pos, int max_string)
  {
    int temp = ceil((double)numbers / 64.);
    uint64_t *bitmap = new uint64_t[temp];
    memcpy(bitmap, bitmap_pos, temp * 8);

    std::vector<long_int> outlier;
    for(int i=0;i<outlier_num;i++){
        mpz_t tmp;
        mpz_init(tmp);
        mpz_import(tmp, sizeof(long_int), -1, sizeof(in[0]), 0, 0, in);
        long_int outlier_(tmp);
        outlier.emplace_back(outlier_);
        in += sizeof(long_int);
        mpz_clear(tmp);
    }
    int left = 0;
    long_int decode = 0;
    uint64_t start = start_byte;
    uint64_t end = 0;
    uint64_t total_bit = l * numbers;
    int writeind = 0;
    end = start + (int)(total_bit / 8);
    long_int *res = out;
    if (total_bit % 8 != 0)
    {
      end++;
    }
    int ones = 0;
    while (start <= end)
    {
      if (writeind >= numbers)
      {
        break;
      }
      while (left >= l)
      {
        long_int tmp = decode & (((long_int)1 << l) - 1);
        long_int tmpval = tmp & (((long_int)1 << (l - 1)) - 1);
        //std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
        if (!(tmp >> (l - 1)))
        {
          tmpval = -tmpval;
        }

        decode = (decode >> l);


        if (((bitmap[writeind / 64] >> (63 - writeind % 64)) & 1))
        {
          ones++;
          //std::cout<<"this is a outlier "<<writeind<<"temp rank is "<<ranktmp<<" value is "<<outlier[ranktmp-1]<<std::endl;
          out[writeind + start_index] = outlier[ones - 1];
          writeind++;
        }
        else
        {
            tmpval += (start_key + (writeind * slope));
            //std::cout<<"this is not a outlier "<<writeind<<" value is "<<tmpval<<std::endl;
            out[writeind + start_index] = tmpval;
            writeind++;
        }

        left -= l;
      }
      decode += ((long_int)in[start] << left);
      start++;
      left += 8;
    }
  }
