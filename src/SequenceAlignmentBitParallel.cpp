/*
 * SequenceAlignmentBitParallel.cpp
 *
 *  Created on: Oct 18, 2013
 *      Author: shu
 */

#include <iostream>
#include "SequenceAlignmentBitParallel.h"

using namespace std;

namespace sequence_alignment_bit_parallel {

SequenceAlignmentBitParallel::SequenceAlignmentBitParallel() :
    match_(0), mismatch_(0), gap_(0) {
}

SequenceAlignmentBitParallel::~SequenceAlignmentBitParallel() {
}

int SequenceAlignmentBitParallel::SetScores(Score match, Score mismatch, Score gap) {
  if (match > 0 && mismatch < 0 && gap < mismatch) {
    match_ = match;
    mismatch_ = mismatch;
    gap_ = gap;

    /*
     Score max_d_h = match - gap;
     Score max_d_v = match - gap;
     Score min_d_h = gap;
     Score min_d_v = gap;

     cout << " \t";
     for (Score d_h = min_d_h; d_h <= max_d_h; ++d_h) {
     cout << d_h << "\t";
     }
     cout << std::endl;

     for (Score d_v = min_d_v; d_v <= max_d_v; ++d_v) {
     cout << d_v << "\t";
     for (Score d_h = min_d_h; d_h <= max_d_h; ++d_h) {
     int formula_id = 3;
     Score new_d_v = gap_;

     Score tmp_d_v = 0;
     tmp_d_v = mismatch_ - d_h;
     if (tmp_d_v > new_d_v) {
     new_d_v = tmp_d_v;
     formula_id = 2;
     }
     tmp_d_v = d_v + gap - d_h;
     if (tmp_d_v > new_d_v) {
     new_d_v = tmp_d_v;
     formula_id = 4;
     }
     cout << formula_id << ")" << "\t";
     }
     cout << endl;
     }
     */
    return 0;
  }
  return 1;
}

void SequenceAlignmentBitParallel::BuildPeq(Char *string, size_t string_length, Word *p_eq) {
  const Word w = 1;
  for (size_t i = 0; i < string_length; ++i) {
    Char c = string[i];
    p_eq[c] |= w << i;
  }
}

SequenceAlignmentBitParallel::Score SequenceAlignmentBitParallel::CalculateAlignmentScore(
    Word *string0_p_eq, Char *string1, size_t string1_length) {
  return 0;
}

} /* namespace sequence_alignment_bit_parallel */
