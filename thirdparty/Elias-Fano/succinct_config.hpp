#pragma once

/* #undef SUCCINCT_USE_LIBCXX */
#ifndef SUCCINCT_USE_LIBCXX
#    define SUCCINCT_USE_LIBCXX 0
#endif	

#define SUCCINCT_USE_INTRINSICS 1
#ifndef SUCCINCT_USE_INTRINSICS
#    define SUCCINCT_USE_INTRINSICS 0
#endif	

/* #undef SUCCINCT_USE_POPCNT */
#ifndef SUCCINCT_USE_POPCNT
#    define SUCCINCT_USE_POPCNT 0
#endif
