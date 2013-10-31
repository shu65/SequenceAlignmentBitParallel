/*
 * SequenceAlignmentBitParallel.h
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

	static const size_t kWordBitLength = 64;

	SequenceAlignmentBitParallel();
	virtual ~SequenceAlignmentBitParallel();

	int SetScores(Score match, Score mismatch, Score gap);
	void BuildPeq(Char *string, size_t string_length, Word *p_eq);
	Score CalculateAlignmentScore(Word *string0_p_eq, size_t string0_length,
			Char *string1, size_t string1_length);

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
			const std::vector<Word>& delta_v_shift,
			const std::vector<Word>& delta_h,
			Word delta_v_not_max_to_boundary_shift);
	Word GetDeltaVShiftInZoneD(Word not_delta_v_min_positions);

	Word GetDeltaHInNotZoneD(Score output_delta_v_id,
			const std::vector<Word>& delta_v_shift,
			const std::vector<Word>& delta_h);

	Word GetDeltaHInZoneD(
			Word not_delta_h_min_positions);

	Score DecodeScore(size_t string0_length, size_t string1_length,
			const std::vector<Word>& delta_h);

	int PrintScoreGapVectors(size_t string0_length,
			const std::vector<Word>& score_gap_vectors);


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
