/*
 * SequenceAlignmentBitParallelTest.cpp
 *
 *   Copyright (c) 2013, Shuji Suzuki
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 *   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */


#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include "../src/SequenceAlignmentBitParallel.h"

TEST(EditDistanceBitParallelTest, CalculateEditDistance) {
  typedef sequence_alignment_bit_parallel::SequenceAlignmentBitParallel TestSequenceAlignmentBitParall;

  TestSequenceAlignmentBitParall::Char alphabet_size = 4;
  TestSequenceAlignmentBitParall::Char str0_values[] = { 0, 1, 3, 2, 1, 0, 0};
  std::vector<TestSequenceAlignmentBitParall::Char> str0(&str0_values[0],
      &str0_values[0] + sizeof(str0_values) / sizeof(str0_values[0]));

  TestSequenceAlignmentBitParall::Char str1_values[] = { 0, 2, 3, 1, 0, 0};
  std::vector<TestSequenceAlignmentBitParall::Char> str1(&str1_values[0],
      &str1_values[0] + sizeof(str1_values) / sizeof(str1_values[0]));
  std::vector<TestSequenceAlignmentBitParall::Word> str0_p_eq(alphabet_size, 0);

  TestSequenceAlignmentBitParall sequence_alignment_bit_parallel;
  sequence_alignment_bit_parallel.SetScores(2, -3, -5);
  sequence_alignment_bit_parallel.BuildPeq(&str0[0], str0.size(), &str0_p_eq[0]);

  TestSequenceAlignmentBitParall::Score score = sequence_alignment_bit_parallel.CalculateAlignmentScore(
      &str0_p_eq[0], &str1[0], str1.size());
  EXPECT_EQ(2, score);
}

