
#ifndef COMMON_H_
#define COMMON_H_

// C headers (sorted)
#include <errno.h>
#include <fcntl.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <sys/mman.h>
#include <sys/resource.h>
#endif
#include <sys/stat.h>
#ifndef _WIN32
#include <sys/time.h>
#endif
#include <sys/types.h>
#include <time.h>
#include <Eigen/Dense>

// C++ headers (sorted)
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <string.h>
#include <stdlib.h>
#include <iso646.h> // mostly for Microsoft compilers
#include <limits.h>
#include <stdint.h> // part of Visual Studio 2010 and better
#include <vector>
#include <cstddef>
#include <getopt.h>
#include <assert.h>
#include <snappy.h>


#ifdef _MSC_VER
#include <iso646.h>
#include <stdint.h>
#include <intrin.h>

#define __attribute__(n)
#define __restrict__ __restrict

#endif

#endif /* COMMON_H_ */
