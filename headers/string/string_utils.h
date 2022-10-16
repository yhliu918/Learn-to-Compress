#pragma once
#pragma GCC diagnostic ignored "-Wshift-count-overflow"
#include <string.h>
#include "../common.h"
#include <boost/multiprecision/gmp.hpp>


template <typename T>
inline uint32_t bits_T(T v)
{
  uint32_t r(0);
  int length = sizeof(T) * 8;
  if (length > 255 && v >= ((T)1 << (uint8_t)255))
  {
    v >>= 255;
    //v = v/2;
    r += 256;
  }
  if (length > 127 && v >= ((T)1 << (uint8_t)127))
  {
    v >>= 128;
    r += 128;
  }
  if (length > 63 && v >= ((T)1 << (uint8_t)63))
  {
    v >>= 64;
    r += 64;
  }
  if (length > 31 && v >= ((T)1 << (uint8_t)31))
  {
    v >>= 32;
    r += 32;
  }
  if (length > 15 && v >= ((T)1 << (uint8_t)15))
  {
    v >>= 16;
    r += 16;
  }
  if (length > 7 && v >= ((T)1 << (uint8_t)7))
  {
    v >>= 8;
    r += 8;
  }
  if (length > 3 && v >= ((T)1 << (uint8_t)3))
  {
    v >>= 4;
    r += 4;
  }
  if (length > 1 && v >= ((T)1 << (uint8_t)1))
  {
    v >>= 2;
    r += 2;
  }
  if (v >= (T)1)
  {
    r += 1;
  }
  
  return r;
}



long_int convertToLongInt(std::string letter)
{
  std::string s2 = "";
  long_int result = 0;
  int shift = (letter.size() - 1) * 8;

  for (auto i = 0; i < letter.size(); i++)
  {
    char x = letter.at(i);
    result += long_int(x) << shift;
    // std::cout<<x<<" "<<long_int(x)<<" "<<(long_int(x)<<shift)<<" "<<result<<std::endl;
    shift -= 8;
  }
  return result;
}

template <typename T>
inline T convertToASCII(std::string& letter)
{
  std::string s2 = "";
  T result = 0;
  int shift = (letter.size() - 1) * 8;

  for (auto x:letter)
  {
    result += T(x) << (uint8_t)shift;
    shift -= 8;
  }
  return result;
}


// *******************  SUBSET  *******************
long_int convertToLongInt_subset(int min_ascii, int max_ascii, std::string letter)
{
  int set_size = max_ascii - min_ascii + 1;
  std::string s2 = "";
  long_int result = 0;
  for (auto i = 0; i < letter.size(); i++)
  {
    char x = letter.at(i);

    if(long_int(x)-min_ascii<long_int(set_size-1)){
      result *= set_size;
      result += (long_int(x)-min_ascii);
    }
    else{
      result *= set_size;
      result += (long_int(set_size-1));
    }
    
    // std::cout<<x<<" "<<long_int(x)<<" "<<(long_int(x)<<shift)<<" "<<result<<std::endl;
  }
  return result;
}


template <typename T>
inline T convertToASCII_subset(int min_ascii, int max_ascii,std::string& letter)
{
  int set_size = max_ascii - min_ascii + 1;
  std::string s2 = "";
  T result = 0;

  for (auto x:letter)
  {
    result *=set_size;
    result += std::min(T(x)-min_ascii, T(set_size-1));
  }
  return result;
}


template <typename T>
inline std::string convertToString_subset(int min_ascii, int max_ascii,T *record, int str_len)
{
  // std::cout<<sizeof(record)<<std::endl;
  int set_size = max_ascii - min_ascii + 1;
  char* res = new char[str_len];
  int i = str_len;
  while((*record)>0){
    i--;
    T result = (*record)/set_size;
		uint32_t m = (*record) - result*set_size;
    res[i]= (char)(m+min_ascii);
		(*record)=result;

	}
  while(i){
    i--;
    res[i]= (char)(min_ascii);
  }
  std::string val = std::string(res , str_len);
  return val;

}


// *******************  SUBSET_SHIFT  *******************


template <typename T>
inline std::string convertToString_subset_shift(int min_ascii, int max_ascii,T *record, int str_len)
{
  int set_size = max_ascii - min_ascii + 1;
  int shift_size = bits_T<uint32_t>(set_size) -1;
  char* res = new char[str_len];
  int shift = (str_len-1) * shift_size;
  uint32_t mask = (1<<shift_size)-1;
  int i = str_len;
  while(i){
    i--;
    res[i]= (uint32_t)((*record) & mask) + min_ascii;
    (*record) >>= shift_size;
  }

  std::string val = std::string(res , str_len);
  return val;

}





void convertToASCII_mpz(std::string letter, mpz_t *result)
{
  mpz_t tmp;
  mpz_init(tmp);

  int shift = (letter.size() - 1) * 8;

  for (auto i = 0; i < letter.size(); i++)
  {
    int x = letter.at(i);
    // covert int x to mpz_t
    mpz_set_ui(tmp, x);
    mpz_mul_2exp(tmp, tmp, shift);
    mpz_add(*result, *result, tmp);
    // std::cout<<x<<" "<<*result<<std::endl;
    shift -= 8;
  }

  mpz_clear(tmp);
}


std::string convertToString_long_int(long_int *record)
{
  // std::cout<<sizeof(record)<<std::endl;
  std::string val;
  int k = 0;
  int len = (mpz_sizeinbase(record->backend().data(), 2) + 63) / 64;
  int str_len = len * 8;
  val.resize(str_len);
  char *input = reinterpret_cast<char *>(&(record->backend().data()->_mp_d[0]));
  int i = 0;
  for (int j = 0; j < 8; ++j)
  {
    if (input[str_len - 1 - j] != 0)
      val[i++] = input[str_len - 1 - j];
  }
  for (int j = 8; j < str_len; ++j)
  {
    val[i++] = input[str_len - 1 - j];
  }
  val.resize(i);
  return val;
}


template <typename T>
inline std::string convertToString(T *record, int str_len)
{
  // std::cout<<sizeof(record)<<std::endl;
  const char * res = reinterpret_cast<const char *>(record);
  std::string val = std::string(res , str_len);
  reverse(val.begin(), val.end());
  //std::string ret(val.rbegin(),val.rend());
  return val;

}



std::string convertToString_mpz(mpz_t *record)
{
  // std::cout<<sizeof(record)<<std::endl;
  std::string val = "";
  int len = (mpz_sizeinbase(*record, 2) + 63) / 64;
  while (len)
  {
    // copy the value of record to uint64 result
    uint64_t result = 0;
    mpz_export(&result, NULL, -1, sizeof(uint64_t), 0, 0, *record);
    // std::cout<<result<<std::endl;

    std::string s = std::string(reinterpret_cast<char *>(&result));
    s = s.substr(0, 8);
    reverse(s.begin(), s.end());
    val = s + val;
    mpz_div_2exp(*record, *record, 64);
    len--;
  }

  return val;
}

int random(int m)
{
    return rand() % m;
}


uint32_t bits_long(long_int v)
{
  uint32_t r(0);
  if (v >= ((long_int)1 << 255))
  {
    v >>= 256;
    r += 256;
  }
  if (v >= ((long_int)1 << 127))
  {
    v >>= 128;
    r += 128;
  }
  if (v >= ((long_int)1 << 63))
  {
    v >>= 64;
    r += 64;
  }
  if (v >= ((long_int)1 << 31))
  {
    v >>= 32;
    r += 32;
  }
  if (v >= ((long_int)1 << 15))
  {
    v >>= 16;
    r += 16;
  }
  if (v >= ((long_int)1 << 7))
  {
    v >>= 8;
    r += 8;
  }
  if (v >= ((long_int)1 << 3))
  {
    v >>= 4;
    r += 4;
  }
  if (v >= ((long_int)1 << 1))
  {
    v >>= 2;
    r += 2;
  }
  if (v >= (long_int)1)
  {
    r += 1;
  }
  return r;
}

uint32_t bits_mpz(mpz_t *v)
{
  uint32_t r(0);
  mpz_t imm_1, shift1, shift3, shift7, shift15, shift31, shift63, shift127, shift255;
  mpz_init(imm_1);
  mpz_init(shift1);
  mpz_init(shift3);
  mpz_init(shift7);
  mpz_init(shift15);
  mpz_init(shift31);
  mpz_init(shift63);
  mpz_init(shift127);
  mpz_init(shift255);

  mpz_set_ui(imm_1, 1);
  mpz_mul_2exp(shift255, imm_1, 255);
  mpz_mul_2exp(shift127, imm_1, 127);
  mpz_mul_2exp(shift63, imm_1, 63);
  mpz_mul_2exp(shift31, imm_1, 31);
  mpz_mul_2exp(shift15, imm_1, 15);
  mpz_mul_2exp(shift7, imm_1, 7);
  mpz_mul_2exp(shift3, imm_1, 3);
  mpz_mul_2exp(shift1, imm_1, 1);

  if (mpz_cmp(*v, shift255) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 256);
    r += 256;
  }
  if (mpz_cmp(*v, shift127) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 128);
    r += 128;
  }
  if (mpz_cmp(*v, shift63) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 64);
    r += 64;
  }
  if (mpz_cmp(*v, shift31) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 32);
    r += 32;
  }
  if (mpz_cmp(*v, shift15) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 16);
    r += 16;
  }
  if (mpz_cmp(*v, shift7) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 8);
    r += 8;
  }
  if (mpz_cmp(*v, shift3) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 4);
    r += 4;
  }
  if (mpz_cmp(*v, shift1) >= 0)
  {
    mpz_tdiv_q_2exp(*v, *v, 2);
    r += 2;
  }
  if (mpz_cmp(*v, imm_1) >= 0)
  {
    r += 1;
  }
  return r;
}
