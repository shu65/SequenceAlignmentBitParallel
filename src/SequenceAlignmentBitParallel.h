/*
 * SequenceAlignmentBitParallel.h
 *
 *  Created on: Oct 18, 2013
 *      Author: shu
 */

#ifndef SEQUENCEALIGNMENTBITPARALLEL_H_
#define SEQUENCEALIGNMENTBITPARALLEL_H_

#include <cstdlib>
#include <limits.h>
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
  Score CalculateAlignmentScore(Word *string0_p_eq, Char *string1,
      size_t string1_length);

private:
  struct InputPair {
    uint32_t delta_v_id;
    uint32_t delta_h_id;
  };
  static const Word kAllOneWord = UINT_MAX;
  size_t GetScoreGapVectorId(Score score) {
    return score - min_score_gap_;
  }
  void BuildOutputInputLists();
  Word GetDeltaVMaxShift(Word Matches, Word delta_h_min);
  Word GetDeltaVShiftInOtherZoneA(const Score output_delta_v_id,
      const Word remain_delta_h_min,
      const std::vector<Word> &delta_v_shift_relative_matches,
      const std::vector<Word> &delta_h);
  Word GetDeltaVShiftInZoenBC(Score output_delta_v_id,
      const std::vector<Word>& delta_v_shift, const std::vector<Word>& delta_h,
      Word delta_v_not_max_to_boundary_shift);
  Word GetDeltaVShiftInZoneD(const size_t score_gap_min_id, const size_t score_gap_max_id,
      const std::vector<Word>& delta_v_shift);

  Word GetDeltaH(Score output_delta_v_id, const std::vector<Word>& delta_v_shift, const std::vector<Word>& delta_h);
  Score DecodeScore(size_t string1_length, const std::vector<Word>& delta_h);

  int PrintScoreGapVectors(const std::vector<Word>& score_gap_vectors);

  Score match_;
  Score mismatch_;
  Score gap_;
  Score min_score_gap_;
  Score max_score_gap_;
  std::vector<std::vector<InputPair> > delta_v_input_pairs_lists_in_a_;
  std::vector<std::vector<InputPair> > delta_v_input_pairs_lists_in_b_;
  std::vector<size_t> delta_v_input_delta_h_list_in_c_;
  std::vector<std::vector<InputPair> > delta_h_input_pairs_lists_;

};

} /* namespace sequence_alignment_bit_parallel */

#endif /* SEQUENCEALIGNMENTBITPARALLEL_H_ */
