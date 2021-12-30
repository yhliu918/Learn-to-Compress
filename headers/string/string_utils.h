#pragma once
#include <string.h>
#include "../common.h"
#include <boost/multiprecision/gmp.hpp>
long_int convertToASCII(std::string letter)
{
    std::string s2= "";
    for (auto i = 0; i < letter.length(); i++)
    {
        char x = letter.at(i);
        std::string stringx = std::to_string(int(x));
        while (stringx.size()<3){
            stringx = '0'+stringx;
        }
        s2+=stringx;
    }
    while(s2[0]=='0'){
        s2 = s2.substr(1,s2.size());
    }
    long_int result{s2.c_str()};
    return result;
}

std::string convertToString(long_int  record){
    std::string result = record.convert_to<std::string>();
    std::string val = "";
    while(result.size()%3!=0){
      result = '0'+result;
    }
    for(int i=0;i<result.size();i+=3){
        std::string tmp_string = result.substr(i,3);
        char tmp = std::stoi(tmp_string);
        val = val+tmp;
    }
    std::cout<<val<<std::endl;
    return val;
    
}


uint32_t bits_long(long_int v) {
  uint32_t r(0);
  if (v >= ((long_int)1 << 255)) {
    v >>= 256;
    r += 256;
  } 
  if (v >= ((long_int)1 << 127)) {
    v >>= 128;
    r += 128;
  }
  if (v >= ((long_int)1 << 63)) {
    v >>= 64;
    r += 64;
  }
  if (v >= ((long_int)1 << 31)) {
    v >>= 32;
    r += 32;
  }
  if (v >= ((long_int)1 << 15)) {
    v >>= 16;
    r += 16;
  }
  if (v >= ((long_int)1 << 7)) {
    v >>= 8;
    r += 8;
  }
  if (v >= ((long_int)1 << 3)) {
    v >>= 4;
    r += 4;
  }
  if (v >= ((long_int)1 << 1)) {
    v >>= 2;
    r += 2;
  }
  if (v >= (long_int)1) {
    r += 1;
  }
  return r;
}

