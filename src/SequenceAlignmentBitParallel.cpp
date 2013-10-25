/*
 * SequenceAlignmentBitParallel.cpp
 *
 *  Created on: Oct 18, 2013
 *      Author: shu
 */

#include <iostream>
#include <vector>
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
    max_score_gap_ = match_ - gap_;
    min_score_gap_ = gap_;

    cout << " \t";
    for (Score d_h = min_score_gap_; d_h <= max_score_gap_; ++d_h) {
      cout << d_h << "\t";
    }
    cout << std::endl;
    for (Score v = min_score_gap_; v <= max_score_gap_; ++v) {
      cout << v << "\t";
      for (Score h = min_score_gap_; h <= max_score_gap_; ++h) {
        int formula_id = 3;
        Score new_d_v = gap_;

        Score tmp_d_v = 0;
        tmp_d_v = mismatch_ - h;
        if (tmp_d_v > new_d_v) {
          new_d_v = tmp_d_v;
          formula_id = 2;
        }
        tmp_d_v = v + gap - h;
        if (tmp_d_v > new_d_v) {
          new_d_v = tmp_d_v;
          formula_id = 4;
        }
        cout << new_d_v + h - v << "(" << new_d_v << ")" << "\t";

      }
      cout << endl;
    }

    /*
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

  vector<Word> delta_v(max_score_gap_ - min_score_gap_, 0);
  vector<Word> delta_h(max_score_gap_ - min_score_gap_, 0);
  delta_h[GetDeltaVectorId(min_score_gap_)] = kAllOneWord;
  for (size_t i = 0; i < string1_length; ++i) {
    delta_v[GetDeltaVectorId(min_score_gap_)] = 1;
    Word matches = string0_p_eq[string1[i]];
    Word delta_v_max_shift = GetDeltaVMaxShift(matches, delta_h[GetDeltaVectorId(min_score_gap_)]);
    Word remain_delta_h_min = delta_h[GetDeltaVectorId(min_score_gap_)] ^ (delta_v_max_shift >> 1);

    Word INITpos6s = delta_h[GetDeltaVectorId(-4)] & (delta_v_max_shift | matches);
    Word DVpos6shift = ((INITpos6s << 1) + remain_delta_h_min) ^ remain_delta_h_min;
    Word DVpos6shiftNotMatch = DVpos6shift & !matches;

    Word INITpos5s = (delta_h[GetDeltaVectorId(-3)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(-4)] & DVpos6shiftNotMatch);
    Word DVpos5shift = ((INITpos6s << 1) + remain_delta_h_min) ^ remain_delta_h_min;
    Word DVpos5shiftNotMatch = DVpos5shift & !matches;

    Word INITpos4s = (delta_h[GetDeltaVectorId(-2)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(-3)] & DVpos6shiftNotMatch)
        | (delta_h[GetDeltaVectorId(-4)] & DVpos5shiftNotMatch);
    Word DVpos4shift = ((INITpos4s << 1) + remain_delta_h_min) ^ remain_delta_h_min;
    Word DVpos4shiftNotMatch = DVpos4shift & !matches;

    Word INITpos3s = (delta_h[GetDeltaVectorId(-1)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(-2)] & DVpos6shiftNotMatch)
        | (delta_h[GetDeltaVectorId(-3)] & DVpos5shiftNotMatch)
        | (delta_h[GetDeltaVectorId(-4)] & DVpos4shiftNotMatch);
    Word DVpos3shift = ((INITpos3s << 1) + remain_delta_h_min) ^ remain_delta_h_min;
    Word DVpos3shiftNotMatch = DVpos4shift & !matches;

    Word DVnot7to3shiftorMatch = ~((delta_v_max_shift | matches) | DVpos6shift | DVpos5shift
        | DVpos4shift | DVpos3shift);
    Word DVpos2shift = ((delta_h[GetDeltaVectorId(0)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(-1)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(-2)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(-3)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(-4)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(-5)] & DVnot7to3shiftorMatch)) << 1;

    Word DVpos1shift = ((delta_h[GetDeltaVectorId(1)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(0)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(-1)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(-2)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(-3)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(-4)] & DVnot7to3shiftorMatch)) << 1;

    Word DV0shift = ((delta_h[GetDeltaVectorId(2)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(1)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(0)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(-1)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(-2)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(-3)] & DVnot7to3shiftorMatch)) << 1;

    Word DVneg1shift = ((delta_h[GetDeltaVectorId(3)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(2)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(1)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(0)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(-1)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(-2)] & DVnot7to3shiftorMatch)) << 1;

    Word DVneg2shift = ((delta_h[GetDeltaVectorId(4)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(3)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(2)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(1)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(0)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(-1)] & DVnot7to3shiftorMatch)) << 1;

    Word DVneg3shift = ((delta_h[GetDeltaVectorId(5)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(4)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(3)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(2)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(1)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(0)] & DVnot7to3shiftorMatch)) << 1;

    Word DVneg4shift = ((delta_h[GetDeltaVectorId(6)] & (delta_v_max_shift | matches))
        | (delta_h[GetDeltaVectorId(5)] & DVpos6shift)
        | (delta_h[GetDeltaVectorId(4)] & DVpos5shift)
        | (delta_h[GetDeltaVectorId(3)] & DVpos4shift)
        | (delta_h[GetDeltaVectorId(2)] & DVpos3shift)
        | (delta_h[GetDeltaVectorId(1)] & DVnot7to3shiftorMatch)) << 1;

    Word DVneg5shift = kAllOneWord
        ^ (delta_v_max_shift | DVpos6shift | DVpos5shift | DVpos4shift | DVpos3shift | DVpos2shift
            | DVpos1shift | DV0shift | DVneg1shift | DVneg2shift | DVneg3shift | DVneg4shift);

  }
  return 0;
}

SequenceAlignmentBitParallel::Word SequenceAlignmentBitParallel::GetDeltaVMaxShift(Word matches,
    Word delta_h_min) {
  Word init_positions = delta_h_min & matches;
  return ((init_positions + delta_h_min) ^ delta_h_min) ^ init_positions;
}

} /* namespace sequence_alignment_bit_parallel */
