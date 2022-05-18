#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <inttypes.h>
#include "common.h"

#define P10_UINT64 10000000000000000000ULL /* 19 zeroes */
#define E10_UINT64 19
#define STRINGIZER(x) #x
#define TO_STRING(x) STRINGIZER(x)

int print_u128_u(uint128_t u128)
{
  int rc;
  if (u128 > UINT64_MAX)
  {
    uint128_t leading = u128 / P10_UINT64;
    uint64_t trailing = u128 % P10_UINT64;
    rc = print_u128_u(leading);
    rc += printf("%." TO_STRING(E10_UINT64) PRIu64, trailing);
  }
  else
  {
    uint64_t u64 = u128;
    rc = printf("%" PRIu64, u64);
  }
  return rc;
}