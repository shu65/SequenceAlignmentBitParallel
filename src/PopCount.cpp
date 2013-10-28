/*
 * PopCount.cpp
 *
 *  Created on: 2013/10/28
 *      Author: shuji
 */

#include <stdint.h>

namespace pop_count {
   uint32_t PopCount64(uint64_t x) {
    x = (x & 0x5555555555555555ULL) + ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
    x = (x & 0x3333333333333333ULL) + ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
    x *= 0x0101010101010101ULL;
    return (uint32_t)x;
  }
}

