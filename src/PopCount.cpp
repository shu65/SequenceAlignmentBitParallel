/*
 * PopCount.cpp
 *
 *  Created on: 2013/10/28
 *      Author: shuji
 */

#include <stdint.h>

namespace pop_count {
uint64_t PopCount64(uint64_t x) {
  x = ((x & 0xaaaaaaaaaaaaaaaaUL) >> 1) + (x & 0x5555555555555555UL);
  x = ((x & 0xccccccccccccccccUL) >> 2) + (x & 0x3333333333333333UL);
  x = ((x & 0xf0f0f0f0f0f0f0f0UL) >> 4) + (x & 0x0f0f0f0f0f0f0f0fUL);
  x = ((x & 0xff00ff00ff00ff00UL) >> 8) + (x & 0x00ff00ff00ff00ffUL);
  x = ((x & 0xffff0000ffff0000UL) >> 16) + (x & 0x0000ffff0000ffffUL);
  x = ((x & 0xffffffff00000000UL) >> 32) + (x & 0x00000000ffffffffUL);
  return x;
}
}

