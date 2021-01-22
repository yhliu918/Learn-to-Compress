/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 */

#ifndef TIME
#define TIME

#include "common.h"

namespace Codecset {

double getNow() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec / 1000000.0;
}
    

} // namespace FastPFor

#endif
