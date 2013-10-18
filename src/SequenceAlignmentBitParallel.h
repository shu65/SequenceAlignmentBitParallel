/*
 * SequenceAlignmentBitParallel.h
 *
 *  Created on: Oct 18, 2013
 *      Author: shu
 */

#ifndef SEQUENCEALIGNMENTBITPARALLEL_H_
#define SEQUENCEALIGNMENTBITPARALLEL_H_

#include <cstdlib>
#include <stdint.h>

namespace sequence_alignment_bit_parallel {

class SequenceAlignmentBitParallel {
public:
  typedef int32_t Score;
  typedef uint8_t Char;
  typedef uint64_t Word;

  SequenceAlignmentBitParallel();
  virtual ~SequenceAlignmentBitParallel();

  int SetScores(Score match, Score mismatch, Score gap);
  void BuildPeq(Char *string, size_t string_length, Word *p_eq);
  Score CalculateAlignmentScore(Word *string0_p_eq, Char *string1, size_t string1_length);

private:
  Score match_;
  Score mismatch_;
  Score gap_;
};

} /* namespace sequence_alignment_bit_parallel */

#endif /* SEQUENCEALIGNMENTBITPARALLEL_H_ */