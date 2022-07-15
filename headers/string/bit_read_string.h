
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
#include "print128.h"

void read_all_bit_fix_string(uint8_t *in, int start_byte, int start_index, int numbers, int l, long_int &slope, long_int &start_key, std::vector<std::string> &result)
{

  int left = 0;
  // long_int decode = 0;
  mpz_t decode_mpz;
  mpz_init(decode_mpz);

  uint64_t start = start_byte;
  uint64_t end = 0;
  uint64_t total_bit = l * numbers;
  int writeind = 0;
  end = start + (int)(total_bit / 8);

  mpz_t tmp_mpz;
  mpz_init(tmp_mpz);
  mpz_t tmp2_mpz;
  mpz_init(tmp2_mpz);
  mpz_t tmpval_mpz;
  mpz_init(tmpval_mpz);
  mpz_t in_tmp;
  mpz_init(in_tmp);

  if (total_bit % 8 != 0)
  {
    end++;
  }
  while (start <= end)
  {

    while (left >= l)
    {
      // long_int tmp = decode & (((long_int)1 << l) - 1);
      mpz_tdiv_r_2exp(tmp_mpz, decode_mpz, l);

      // long_int tmpval = tmp & (((long_int)1 << (l - 1)) - 1);
      mpz_tdiv_r_2exp(tmpval_mpz, tmp_mpz, l - 1);

      if (!mpz_tstbit(tmp_mpz, l - 1))
      {
        // tmpval = -tmpval;
        mpz_neg(tmpval_mpz, tmpval_mpz);
      }
      // decode = (decode >> l);
      mpz_tdiv_q_2exp(decode_mpz, decode_mpz, l);

      // tmpval += (start_key + writeind * slope);
      mpz_add(tmpval_mpz, tmpval_mpz, start_key.backend().data());
      mpz_addmul_ui(tmpval_mpz, slope.backend().data(), writeind);

      // double start_double = Codecset::getNow();
      long_int tmpval(tmpval_mpz);
      result.emplace_back(convertToString(&tmpval,l));
      // double end_double = Codecset::getNow();
      // totaltime += end_double - start_double;

      writeind++;
      left -= l;
      if (left == 0)
      {
        // decode = 0;
        mpz_set_ui(decode_mpz, 0);
      }
    }
    // decode += ((long_int)in[start] << left);

    mpz_set_ui(in_tmp, in[start]);
    mpz_mul_2exp(tmp2_mpz, in_tmp, left);
    mpz_add(decode_mpz, decode_mpz, tmp2_mpz);

    start++;
    left += 8;
  }

  mpz_clear(tmp_mpz);
  mpz_clear(tmp2_mpz);
  mpz_clear(in_tmp);
  mpz_clear(tmpval_mpz);
  mpz_clear(decode_mpz);
}

void read_all_bit_fix_string_mpz(uint8_t *in, int start_byte, int start_index, int numbers, int l, mpz_t *slope, mpz_t *start_key, std::vector<mpz_t *> out)
{

  int left = 0;
  mpz_t decode;
  mpz_init(decode);
  uint64_t start = start_byte;
  uint64_t end = 0;
  uint64_t total_bit = l * numbers;
  int writeind = 0;
  end = start + (int)(total_bit / 8);

  if (total_bit % 8 != 0)
  {
    end++;
  }
  while (start <= end)
  {

    while (left >= l)
    {
      if (writeind >= out.size())
      {
        break;
      }
      mpz_t tmp;
      mpz_init(tmp);
      mpz_t tmpval;
      mpz_init(tmpval);
      // tmp = decode & ((1 << l) - 1);

      // tmpval = tmp & ((1 << (l - 1)) - 1);
      // decode and with (1 right shift l-1)
      mpz_t imm1;
      mpz_init(imm1);
      mpz_set_ui(imm1, 1);

      mpz_t imml;
      mpz_init(imml);
      mpz_mul_2exp(imml, imm1, l);
      mpz_sub_ui(imml, imml, 1);

      mpz_t imml_1;
      mpz_init(imml_1);
      mpz_mul_2exp(imml_1, imm1, l - 1);
      mpz_sub_ui(imml_1, imml_1, 1);

      mpz_and(tmp, decode, imml);
      mpz_and(tmpval, tmp, imml_1);

      // if (!(((tmp >> (l - 1)) & (long_int)1)))
      if (mpz_tstbit(tmp, l - 1) == 0)
      {
        // tmpval = -tmpval;
        mpz_neg(tmpval, tmpval);
      }

      // decode right shift l, decode = decode >> l
      mpz_fdiv_q_2exp(decode, decode, l);
      mpz_add(tmpval, tmpval, *start_key);
      mpz_addmul_ui(tmpval, *slope, writeind);
      // std::cout<<"index "<<writeind<<" tmpval "<<tmpval<<std::endl;
      mpz_set(*out[writeind], tmpval);
      writeind++;
      left -= l;
      if (left == 0)
      {
        mpz_set_ui(decode, 0);
      }

      mpz_clear(tmp);
      mpz_clear(tmpval);
      mpz_clear(imm1);
      mpz_clear(imml);
      mpz_clear(imml_1);

      // std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
    }
    // decode += (in[start] << left)
    mpz_t in_start;
    mpz_init(in_start);
    // set in_start to in[start]
    mpz_import(in_start, 1, 1, sizeof(uint8_t), 0, 0, in + start);
    // mpz_set_ui(in_start, in[start]);
    mpz_mul_2exp(in_start, in_start, left);
    mpz_add(decode, decode, in_start);

    mpz_clear(in_start);

    start++;
    left += 8;
  }
}

template <typename T>
void read_all_fix_string(uint8_t *in, int start_byte, int start_index, int numbers, int l, T slope, T start_key, std::vector<std::string> &result)
{
  uint64_t left = 0;
  T decode = 0;
  uint64_t start = start_byte;
  uint64_t end = 0;
  uint64_t total_bit = l * numbers;
  uint64_t writeind = 0;
  T mask0 = ((T)1 << l) - 1;
  T mask1 = ((T)1<< (l - 1)) - 1;
  end = start + (uint64_t)(total_bit / 8);

  if (total_bit % 8 != 0)
  {
    end++;
  }

  while (start <= end)
  {
    /*
    if(writeind>= numbers){
      break;
    }
    */
    while (left >= l)
    {
      T tmp = decode & mask0;
      T tmpval = tmp & mask1;
      T record_val = 0;
      bool sign = (tmp >> (l - 1)) & 1;

      decode = (decode >> l);
      if(sign){
          record_val = tmpval + (T)(start_key + writeind * slope);
      }
      else{
          record_val = (T)(start_key + writeind * slope) - tmpval;
      }
      
      result[writeind] = convertToString<T>(&record_val, l);
      writeind++;
      left -= l;
      if (left == 0)
      {
        decode = 0;
      }
      // std::cout<<"decode "<<decode<<"left"<<left<<std::endl;
    }
    decode += ((T)in[start] << left);

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
  for (int i = 0; i < outlier_num; i++)
  {
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
      // std::cout<<"decode: "<<decode<<"tmp: "<<tmp<<" l: "<<l<<std::endl;
      if (!(tmp >> (l - 1)))
      {
        tmpval = -tmpval;
      }

      decode = (decode >> l);

      if (((bitmap[writeind / 64] >> (63 - writeind % 64)) & 1))
      {
        ones++;
        // std::cout<<"this is a outlier "<<writeind<<"temp rank is "<<ranktmp<<" value is "<<outlier[ranktmp-1]<<std::endl;
        out[writeind + start_index] = outlier[ones - 1];
        writeind++;
      }
      else
      {
        tmpval += (start_key + (writeind * slope));
        // std::cout<<"this is not a outlier "<<writeind<<" value is "<<tmpval<<std::endl;
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

std::string read_bit_fix_string_(uint8_t *in, uint32_t l, int to_find, long_int theta1, long_int theta0)
{
  uint64_t find_bit = to_find * l;
  uint64_t start_byte = find_bit / 8;
  uint64_t start_bit = find_bit % 8;
  uint64_t occupy = start_bit;

  int total = 0;
  mpz_t in_start;
  mpz_init(in_start);
  mpz_t decode_mpz;
  mpz_init(decode_mpz);
  mpz_t tmp_mpz;
  mpz_init(tmp_mpz);
  mpz_t tmpval_mpz;
  mpz_init(tmpval_mpz);

  while (total < l)
  {
    // decode += ((long_int)(val >> occupy) << total);
    mpz_set_ui(in_start, in[start_byte] >> occupy);
    mpz_mul_2exp(in_start, in_start, total);
    mpz_add(decode_mpz, decode_mpz, in_start);

    total += (8 - occupy);
    occupy = 0;
    start_byte++;
  }
  // long_int tmp = decode & (((long_int)1 << l) - 1);
  mpz_tdiv_r_2exp(tmp_mpz, decode_mpz, l);

  // long_int tmpval = tmp & (((long_int)1 << (l - 1)) - 1);
  mpz_tdiv_r_2exp(tmpval_mpz, tmp_mpz, l - 1);
  if (!mpz_tstbit(tmp_mpz, l - 1))
  {
    // tmpval = -tmpval;
    mpz_neg(tmpval_mpz, tmpval_mpz);
  }
  // long_int out = value + (theta0 + to_find * theta1);
  mpz_add(tmpval_mpz, tmpval_mpz, theta0.backend().data());
  mpz_addmul_ui(tmpval_mpz, theta1.backend().data(), to_find);
  long_int out(tmpval_mpz);
  mpz_clear(tmpval_mpz);
  mpz_clear(tmp_mpz);
  mpz_clear(decode_mpz);
  mpz_clear(in_start);

  return convertToString(&out,0);
}

void read_bit_fix_string_long_int(uint8_t *in, uint32_t l, int to_find, long_int theta1, long_int theta0, mpz_t *result)
{
  uint64_t find_bit = to_find * l;
  uint64_t start_byte = find_bit / 8;
  uint64_t start_bit = find_bit % 8;
  uint64_t occupy = start_bit;

  int total = 0;
  mpz_t in_start;
  mpz_init(in_start);
  mpz_t decode_mpz;
  mpz_init(decode_mpz);
  mpz_t tmp_mpz;
  mpz_init(tmp_mpz);
  mpz_t tmpval_mpz;
  mpz_init(tmpval_mpz);

  while (total < l)
  {
    // decode += ((long_int)(val >> occupy) << total);
    mpz_set_ui(in_start, in[start_byte] >> occupy);
    mpz_mul_2exp(in_start, in_start, total);
    mpz_add(decode_mpz, decode_mpz, in_start);

    total += (8 - occupy);
    occupy = 0;
    start_byte++;
  }
  // long_int tmp = decode & (((long_int)1 << l) - 1);
  mpz_tdiv_r_2exp(tmp_mpz, decode_mpz, l);

  // long_int tmpval = tmp & (((long_int)1 << (l - 1)) - 1);
  mpz_tdiv_r_2exp(tmpval_mpz, tmp_mpz, l - 1);
  if (!mpz_tstbit(tmp_mpz, l - 1))
  {
    // tmpval = -tmpval;
    mpz_neg(tmpval_mpz, tmpval_mpz);
  }
  // long_int out = value + (theta0 + to_find * theta1);
  mpz_add(tmpval_mpz, tmpval_mpz, theta0.backend().data());
  mpz_addmul_ui(tmpval_mpz, theta1.backend().data(), to_find);
  // set out to tmpval_mpz
  // mpz_clear(tmpval_mpz);
  mpz_clear(tmp_mpz);
  mpz_clear(decode_mpz);
  mpz_clear(in_start);
  mpz_set(*result, tmpval_mpz);
}

void read_bit_fix_string_long_int_128(uint8_t *in, uint32_t l, int to_find, uint128_t theta1, uint128_t theta0, uint128_t *result)
{
  uint64_t find_bit = to_find * l;
  uint64_t start_byte = find_bit / 8;
  uint64_t start_bit = find_bit % 8;
  uint64_t occupy = start_bit;
  uint128_t decode = 0;
  uint64_t total = 0;
  while (total < l)
  {
    decode += ((in[start_byte] >> occupy) << total);
    total += (8 - occupy);
    occupy = 0;
    start_byte++;
  }

  decode = decode & (((uint128_t)1 << l) - 1);
  bool sign = (decode & ((uint128_t)1 << (l - 1)));
  int128_t value = (decode & (((uint128_t)1 << (l - 1)) - 1));
  uint128_t tmpval = (decode & (((uint128_t)1 << (l - 1)) - 1));
  // std::cout<<"value:"<<std::endl;
  // print_u128_u(tmpval);
  // std::cout<<std::endl;
  if (!sign)
  {
    value = -value;
  }
  // std::cout<<"l: "<<l<<" value: "<<value<<" predict: "<< (long long) (start_key +((float)to_find*slope))<<std::endl;
  *result = value + (uint128_t)(theta0 + to_find * theta1);
}
template <typename T>
void read_bit_fix_string(const uint8_t *in, uint32_t l, uint32_t to_find, T theta1, T theta0, T *result, uint8_t* ori_length)
{
  uint64_t find_bit = to_find * (l+8);
  uint64_t start_byte = find_bit / 8;
  uint8_t start_bit = find_bit % 8;
  uint64_t occupy = start_bit;
  uint64_t total = 0;

  // decode += (((T)(in[start_byte] >> occupy)) << total);
  // total += (8 - occupy);
  // start_byte++;

  // while (total < l+8)
  // {
  //   decode += ((T)(in[start_byte]) << total);
  //   total += 8 ;
  //   start_byte++;
  // }


  
  T decode = 0;
  memcpy(&decode, in+start_byte, sizeof(T));
  decode >>=(uint8_t)start_bit;
  decode &= (((T)1<<(uint8_t)(l+8))-1);
  // T one = 1;
  // one.left_shift((uint8_t)(l+8) ,*result);
  // decode &= (*result - 1);

  *ori_length = (uint8_t)(decode >> (uint8_t)l); // TODO: maybe right shift?
  bool sign = (decode >> (uint8_t)(l - 1)) & 1;
  T value = (decode & (((T)1 << (uint8_t)(l - 1)) - 1));
  T pred = theta0 + theta1 * to_find;
  if (sign)
  {
    *result = value + pred;

  }
  else
  {
    *result = pred - value;

  }
}

