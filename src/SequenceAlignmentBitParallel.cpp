/*
 * SequenceAlignmentBitParallel.cpp
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

#include <iostream>
#include <vector>
#include <assert.h>
#include "pop_count/PopCount.h"
#include "SequenceAlignmentBitParallel.h"

using namespace std;

namespace sequence_alignment_bit_parallel {

SequenceAlignmentBitParallel::SequenceAlignmentBitParallel() :
		match_(0), mismatch_(0), gap_(0), min_score_gap_(gap_), max_score_gap_(
				match_ - gap_) {
}

SequenceAlignmentBitParallel::~SequenceAlignmentBitParallel() {
}

int SequenceAlignmentBitParallel::SetScores(Score match, Score mismatch,
		Score gap) {
	if (match > 0 && mismatch < 0 && gap < mismatch) {
		match_ = match;
		mismatch_ = mismatch;
		gap_ = gap;
		min_score_gap_ = gap_;
		max_score_gap_ = match_ - gap_;
		BuildOutputInputLists();
		return 0;
	}
	return 1;
}

void SequenceAlignmentBitParallel::BuildPeq(Char *string, size_t string_length,
		Word *p_eq) {
	const Word w = 1;
	for (size_t i = 0; i < string_length; ++i) {
		Char c = string[i];
		p_eq[c] |= w << i;
	}
}

SequenceAlignmentBitParallel::Score SequenceAlignmentBitParallel::CalculateAlignmentScore(
		Word *string0_p_eq, size_t string0_length, Char *string1,
		size_t string1_length) {
	const size_t score_gap_min_id = GetScoreGapVectorId(min_score_gap_);
	const size_t score_gap_max_id = GetScoreGapVectorId(max_score_gap_);
	const Score min_score_in_zone_a = mismatch_ - gap_ + 1;
	const size_t score_gap_boundary_id = GetScoreGapVectorId(
			min_score_in_zone_a);
	const size_t score_gap_vector_length = max_score_gap_ - min_score_gap_ + 1;
	vector<Word> delta_v_shift(score_gap_vector_length, 0);
	vector<Word> delta_v_shift_relative_matches(score_gap_vector_length, 0);
	vector<vector<Word> > delta_hs(2);
	for (size_t i = 0; i < delta_hs.size(); ++i) {
		delta_hs[i].resize(score_gap_vector_length, 0);
	}
	vector<Word> delta_h_relative_matches(score_gap_vector_length, 0);
	size_t delta_hs_read_id = 0;
	size_t delta_hs_write_id = 1 - delta_hs_read_id;
	delta_hs[delta_hs_read_id][min_score_gap_] = kAllOneWord;
	for (size_t i = 0; i < string1_length; ++i, delta_hs_read_id =
			delta_hs_write_id, delta_hs_write_id = 1 - delta_hs_write_id) {
		const vector<Word> &delta_h = delta_hs[delta_hs_read_id];
#ifdef DEBUG
		cout << i << ":\t";
		assert(!PrintScoreGapVectors(string0_length, delta_hs[delta_hs_read_id]));
#endif
		Word matches = string0_p_eq[string1[i]];
		Word not_matches = ~matches;

		delta_v_shift[score_gap_max_id] = GetDeltaVMaxShift(matches,
				delta_h[score_gap_min_id]);
		Word remain_delta_h_min = delta_h[score_gap_min_id]
				^ (delta_v_shift[score_gap_max_id] >> 1);
		delta_v_shift_relative_matches[score_gap_max_id] =
				delta_v_shift[score_gap_max_id] | matches;
		Word delta_v_not_max_to_boundary_shift =
				delta_v_shift_relative_matches[score_gap_max_id];

		for (size_t output_score_id = score_gap_max_id - 1;
				output_score_id >= score_gap_boundary_id; --output_score_id) {
			delta_v_shift[output_score_id] = GetDeltaVShiftInOtherZoneA(
					output_score_id, remain_delta_h_min,
					delta_v_shift_relative_matches, delta_h);
			delta_v_shift_relative_matches[output_score_id] =
					delta_v_shift[output_score_id] & not_matches;
			delta_v_not_max_to_boundary_shift |= delta_v_shift[output_score_id];
		}
		delta_v_not_max_to_boundary_shift = ~delta_v_not_max_to_boundary_shift;

		for (size_t output_score_id = score_gap_boundary_id - 1;
				output_score_id > score_gap_min_id; --output_score_id) {
			delta_v_shift[output_score_id] = GetDeltaVShiftInZoenBC(
					output_score_id, delta_v_shift_relative_matches, delta_h,
					delta_v_not_max_to_boundary_shift);
		}
		delta_v_shift[score_gap_min_id] = GetDeltaVShiftInZoneD(
				score_gap_min_id, score_gap_max_id, delta_v_shift);

#ifdef DEBUG
		//cout << "v:\t";
		//assert(!PrintScoreGapVectors(string0_length, delta_v_shift));
#endif
		vector<Word> &new_delta_h = delta_hs[delta_hs_write_id];

		for (size_t i = score_gap_min_id; i < score_gap_max_id; ++i) {
			delta_h_relative_matches[i] = delta_h[i] & not_matches;
		}
		delta_h_relative_matches[score_gap_max_id] = delta_h[score_gap_max_id]
				| matches;

		Word not_delta_h_min_positions = 0;
		for (size_t output_score_id = score_gap_min_id + 1;
				output_score_id <= score_gap_max_id; ++output_score_id) {
			new_delta_h[output_score_id] = GetDeltaH(output_score_id,
					delta_v_shift, delta_h_relative_matches);
			not_delta_h_min_positions |= new_delta_h[output_score_id];
		}
		new_delta_h[score_gap_min_id] = kAllOneWord ^ not_delta_h_min_positions;
	}
#ifdef DEBUG
	cout << string1_length << ":\t";
	assert(!PrintScoreGapVectors(string0_length, delta_hs[delta_hs_read_id]));
#endif
	return DecodeScore(string0_length, string1_length,
			delta_hs[delta_hs_read_id]);
}

void SequenceAlignmentBitParallel::BuildOutputInputLists() {
	delta_v_input_pairs_lists_in_a_.clear();
	delta_v_input_pairs_lists_in_b_.clear();
	delta_v_input_delta_h_list_in_c_.clear();
	delta_h_input_pairs_lists_.clear();

	size_t score_gap_vector_length = max_score_gap_ - min_score_gap_ + 1;
	delta_v_input_pairs_lists_in_a_.resize(score_gap_vector_length);
	delta_v_input_pairs_lists_in_b_.resize(score_gap_vector_length);
	delta_v_input_delta_h_list_in_c_.resize(score_gap_vector_length, 0);
	delta_h_input_pairs_lists_.resize(score_gap_vector_length);
	Score score_gap_min = gap_;
	Score score_gap_boundary = mismatch_ - gap_;
	Score min_score_in_zone_a = mismatch_ - gap_ + 1;
#ifdef DEBUG
	cout << " \t";
	for (Score d_h = min_score_gap_; d_h <= max_score_gap_; ++d_h) {
		cout << d_h << "("<< GetScoreGapVectorId(d_h) << ")\t";
	}
	cout << std::endl;
#endif
	for (Score delta_v = min_score_gap_; delta_v <= max_score_gap_; ++delta_v) {
#ifdef DEBUG
		cout << delta_v << "("<< GetScoreGapVectorId(delta_v) << ")\t";
#endif
		for (Score delta_h = min_score_gap_; delta_h <= max_score_gap_;
				++delta_h) {
			Score output_delta_v = gap_;
			Score tmp_d_v = 0;
			tmp_d_v = mismatch_ - delta_h;
			if (tmp_d_v > output_delta_v) {
				output_delta_v = tmp_d_v;
			}
			tmp_d_v = delta_v + gap_ - delta_h;
			if (tmp_d_v > output_delta_v) {
				output_delta_v = tmp_d_v;
			}
			InputPair input_pair;
			input_pair.delta_v_id = GetScoreGapVectorId(delta_v);
			input_pair.delta_h_id = GetScoreGapVectorId(delta_h);

			if (output_delta_v >= min_score_in_zone_a) {
				if (delta_h != score_gap_min) {
					delta_v_input_pairs_lists_in_a_[GetScoreGapVectorId(
							output_delta_v)].push_back(input_pair);
				}
			} else if (output_delta_v != score_gap_min) {
				if (delta_v > score_gap_boundary) {
					delta_v_input_pairs_lists_in_b_[GetScoreGapVectorId(
							output_delta_v)].push_back(input_pair);
				} else {
					delta_v_input_delta_h_list_in_c_[GetScoreGapVectorId(
							output_delta_v)] = input_pair.delta_h_id;
				}
			}
			Score output_delta_h = output_delta_v + delta_h - delta_v;
			delta_h_input_pairs_lists_[GetScoreGapVectorId(output_delta_h)].push_back(
					input_pair);
#ifdef DEBUG
			cout << output_delta_v << "(" << output_delta_h << ")" << "\t";
#endif

		}
#ifdef DEBUG
		cout << endl;
#endif
	}
}

SequenceAlignmentBitParallel::Word SequenceAlignmentBitParallel::GetDeltaVMaxShift(
		Word matches, Word delta_h_min) {
	Word init_positions = delta_h_min & matches;
	return ((init_positions + delta_h_min) ^ delta_h_min) ^ init_positions;
}

SequenceAlignmentBitParallel::Word SequenceAlignmentBitParallel::GetDeltaVShiftInOtherZoneA(
		const Score output_delta_v_id, const Word remain_delta_h_min,
		const vector<Word> &delta_v_shift_relative_matches,
		const vector<Word> &delta_h) {
	Word init_positions = 0;
	vector<InputPair> &inputs =
			delta_v_input_pairs_lists_in_a_[output_delta_v_id];
	for (vector<InputPair>::const_iterator inputs_it = inputs.begin();
			inputs_it != inputs.end(); ++inputs_it) {
		init_positions |= delta_h[inputs_it->delta_h_id]
				& delta_v_shift_relative_matches[inputs_it->delta_v_id];
	}
	return ((init_positions << 1) + remain_delta_h_min) ^ remain_delta_h_min;
}

SequenceAlignmentBitParallel::Word SequenceAlignmentBitParallel::GetDeltaVShiftInZoenBC(
		Score output_delta_v_id,
		const vector<Word>& delta_v_shift_relative_matches,
		const vector<Word>& delta_h, Word delta_v_not_max_to_boundary_shift) {
	Word init_positions = 0;
	const vector<InputPair> &inputs =
			delta_v_input_pairs_lists_in_b_[output_delta_v_id];
	for (vector<InputPair>::const_iterator inputs_it = inputs.begin();
			inputs_it != inputs.end(); ++inputs_it) {
		init_positions |= delta_h[inputs_it->delta_h_id]
				& delta_v_shift_relative_matches[inputs_it->delta_v_id];
	}
	init_positions |=
			delta_h[delta_v_input_delta_h_list_in_c_[output_delta_v_id]]
					& delta_v_not_max_to_boundary_shift;
	return init_positions << 1;
}

SequenceAlignmentBitParallel::Word SequenceAlignmentBitParallel::GetDeltaVShiftInZoneD(
		const size_t score_gap_min_id, const size_t score_gap_max_id,
		const vector<Word>& delta_v_shift) {
	Word not_min_positions = 0;
	for (size_t i = score_gap_max_id; i > score_gap_min_id; --i) {
		not_min_positions |= delta_v_shift[i];
	}
	return kAllOneWord ^ not_min_positions;
}

SequenceAlignmentBitParallel::Word SequenceAlignmentBitParallel::GetDeltaH(
		Score output_delta_h_id, const std::vector<Word>& delta_v,
		const std::vector<Word>& delta_h_relative_matches) {
	Word positions = 0;
	vector<InputPair> &inputs = delta_h_input_pairs_lists_[output_delta_h_id];
	for (vector<InputPair>::const_iterator inputs_it = inputs.begin();
			inputs_it != inputs.end(); ++inputs_it) {
		positions |= delta_h_relative_matches[inputs_it->delta_h_id]
				& delta_v[inputs_it->delta_v_id];
	}
	return positions;
}

SequenceAlignmentBitParallel::Score SequenceAlignmentBitParallel::DecodeScore(
		size_t string0_length, size_t string1_length,
		const std::vector<Word>& delta_h) {
	Score sum_score = string1_length * gap_;
	Score current_score = min_score_gap_;
	Word mask = kAllOneWord;
	if (kWordBitLength > string0_length) {
		mask = (1 << string0_length) - 1;
	}
	for (std::vector<Word>::const_iterator delta_h_it = delta_h.begin();
			delta_h_it != delta_h.end(); ++delta_h_it, ++current_score) {
		sum_score += pop_count::PopCount64(*delta_h_it & mask) * current_score;
	}
	return sum_score;
}

int SequenceAlignmentBitParallel::PrintScoreGapVectors(size_t string0_length,
		const std::vector<Word>& score_gap_vectors) {
	Word setted_positions = 0;
	vector<Score> decoded_scores(kWordBitLength, 0);
	Score current_score = min_score_gap_;
	Word mask = kAllOneWord;
	if (kWordBitLength > string0_length) {
		mask = (1 << string0_length) - 1;
	}
	for (std::vector<Word>::const_iterator score_gap_vectors_it =
			score_gap_vectors.begin();
			score_gap_vectors_it != score_gap_vectors.end();
			++score_gap_vectors_it, ++current_score) {
		if ((setted_positions & *score_gap_vectors_it) != 0) {
			return 1;
		}
		setted_positions |= *score_gap_vectors_it & mask;
		Word tmp_vector = *score_gap_vectors_it;
		for (size_t i = 0; i < kWordBitLength; ++i, tmp_vector >>= 1) {
			if (tmp_vector & 1) {
				decoded_scores[i] = current_score;
			}
		}
	}

	if (setted_positions != mask) {
		return 1;
	}

	for (size_t i = 0; i < kWordBitLength; ++i) {
		cout << decoded_scores[i] << "\t";
	}
	cout << endl;

	return 0;
}

}
/* namespace sequence_alignment_bit_parallel */
