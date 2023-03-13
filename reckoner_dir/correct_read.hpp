/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#ifndef CORRECT_READ_HPP
#define CORRECT_READ_HPP

#include "Log.h"
#include "define.hpp"
#include "parse_args.hpp"
#include "Log.h"
#include "QueryKMCDb.h"
#include <algorithm>
#include <cassert>
#include <vector>
#include <fstream>



 //----------------------------------------------------------------------
 // C_candidate_path
 //----------------------------------------------------------------------


class Single_mod {
public:
    std::size_t pos;
    char modification;

    bool insertion;
    bool deletion;

    Single_mod() : insertion(false), deletion(false) {}
    Single_mod(std::size_t _pos) : pos(_pos), modification('\0'), insertion(true), deletion(false) {}
    Single_mod(std::size_t _pos, char _modification, bool _deletion = false) : pos(_pos), modification(_modification), insertion(false), deletion(_deletion) {}
};


class C_candidate_path {
public:
    // modified positions and modification results
    std::vector<Single_mod> modifications;

private:
    int read_length_change;
    std::size_t num_insertions;

    double nucleotides_probability;

    // the values are saved to re-rate after choosing extending k-mer quality
    double saved_kmers_quality;
    double saved_covering_kmers_weight;

    double path_rate;

    // in case of attaching a path it stores number of the attached path modifications
    std::size_t beginning_attached_path_modifications;

    void rate_internal(double kmers_quality, double covering_kmers_weight);

public:
    int get_read_length_change() const { return read_length_change; }
    double get_path_rate() const { return path_rate; }

    void add_substitution(std::size_t pos, char symbol, double symbol_probability);

    void add_insertion(std::size_t pos);

    void add_deletion(std::size_t pos, char symbol);

    void rate(double kmers_quality, std::size_t num_covering_kmers);
    
    void attach_extending_rate(double max_extended_kmer_quality);
    void attach_modifications_at_beginning(const C_candidate_path& beginning_path);

    std::size_t get_beginning_attached_path_modifications() const { return beginning_attached_path_modifications; }

    // constructors
    C_candidate_path() :
        read_length_change(0),
        num_insertions(0),
        nucleotides_probability(1.0),
        saved_kmers_quality(0.0),
        saved_covering_kmers_weight(0.0),
        path_rate(0.0),
        beginning_attached_path_modifications(0) {
    };
};



class C_modification_with_quality {
public:
    float quality;
    char modification; // new symbol in case of substitution or deletion

    bool potential_indel;

    bool insertion;
    bool deletion; // is deletion before this position (when we move towards 3' or in the first k-mer) or after this position (when move towards 5')

    C_modification_with_quality() : potential_indel(false), insertion(false), deletion(false) {}
};



//----------------------------------------------------------------------
// Type wrapping bool variable to prevent std::vector<bool> specialization.
//----------------------------------------------------------------------
template<typename Type>
class TypeWrapper {
public:
    TypeWrapper() : data(Type()) {}
    TypeWrapper(const Type& _data) : data(_data) {}
    operator Type() { return data; }

private:
    Type data;
};



template<bool CORRECT_INDEL>
class C_correct_read {
public:
    // variables
    std::size_t read_length;
    std::size_t quality_score_offset;
    std::size_t kmer_length;
    std::size_t long_kmer_length;
    std::size_t max_extension;
    bool accept_filtered_with_long_kmers;

    // constructors

    C_correct_read(std::size_t _quality_score_offset, std::size_t _kmer_length, std::size_t _long_kmer_length, std::size_t _max_extension,
        const std::string& _read_file_name, CKMCFile& _kmc_file, CKMCFile* _kmc_long_file, bool _accept_filtered_with_long_kmers) :
        read_length(0),
        quality_score_offset(_quality_score_offset),
        kmer_length(_kmer_length),
        long_kmer_length(_long_kmer_length),
        max_extension(_max_extension),
        accept_filtered_with_long_kmers(_accept_filtered_with_long_kmers),
        num_corrected_reads(0),
        num_corrected_errors_step1_1(0),
        num_corrected_errors_step1_2(0),
        num_corrected_errors_step1_3(0),
        num_corrected_errors_step2_1(0),
        num_corrected_errors_step2_2(0),
        num_corrected_substs(0),
        num_corrected_ins(0),
        num_corrected_dels(0),
        num_corrected_pairs(0),
        num_first_kmer_corrections(0),
        num_first_kmer_successes(0),
        num_single_deletions(0),
        num_checked_with_long_kmer(0),
        num_filtered_with_long_kmer(0),
        num_instantly_accepted_with_long_kmer(0),
        num_oriented_reads(0),
        c_err(std::cerr),
        query(_kmc_file, _kmc_long_file, static_cast<uint32>(kmer_length), static_cast<uint32>(long_kmer_length)),
        currently_proposed_indels(0) {
        is_low_LUT.resize(MAX_SCORE + 1);
        for (int i = 0; i < static_cast<int>(quality_score_offset) + QS_CUTOFF; ++i) {
            is_low_LUT[i] = true;
        }
        for (int i = static_cast<int>(quality_score_offset) + QS_CUTOFF; i <= MAX_SCORE; ++i) {
            is_low_LUT[i] = false;
        }

        quality_to_probability_LUT.resize(MAX_SCORE + 1);
        for (char asciiValue = 0; asciiValue <= MAX_SCORE; ++asciiValue) {
            char quality = asciiValue - quality_score_offset;
            if (quality < 0) {
                quality_to_probability_LUT[asciiValue] = 1.0;
            }
            else {
                quality_to_probability_LUT[asciiValue] = std::pow(10.0, -quality / 10.0);
            }
        }
    }

    void correct_errors_in_a_read();
    void prepare_corrector();

    std::string header;
    std::string sequence_modified;
    std::string connector;
    std::string quality_score;

    std::size_t num_corrected_reads;

    std::size_t num_corrected_errors_step1_1;
    std::size_t num_corrected_errors_step1_2;
    std::size_t num_corrected_errors_step1_3;
    std::size_t num_corrected_errors_step2_1;
    std::size_t num_corrected_errors_step2_2;

    std::size_t num_corrected_substs;
    std::size_t num_corrected_ins;
    std::size_t num_corrected_dels;
    std::size_t num_corrected_pairs; // pairs of insertion and deletion in a single region
    std::size_t num_first_kmer_corrections;
    std::size_t num_first_kmer_successes;
    std::size_t num_single_deletions;
    std::size_t num_checked_with_long_kmer;
    std::size_t num_filtered_with_long_kmer;
    std::size_t num_instantly_accepted_with_long_kmer;

    std::size_t num_oriented_reads;

private:
    C_log c_err;

    QueryKMCDb query;

    // sequence of changes in read introduced by recursive functions
    std::vector<C_modification_with_quality> modifications_sequence_stack;
    // was a potential indel detected while the current path extending?
    int currently_proposed_indels;
    int max_current_indels = MAX_INDELS;
    int current_region_index_start = 0;

    std::vector<TypeWrapper<bool>> is_low_LUT;
    std::vector<double> quality_to_probability_LUT;

    // functions
    void correct_errors_between_solid_regions(const std::size_t index_start, const std::size_t index_end, std::vector<C_candidate_path>& candidate_path_vector_out);
    void correct_errors_5_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out);
    void correct_errors_3_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out);
    void correct_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector);
    void extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);
    void extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);
    void extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);
    bool perform_extend_out_left(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector);
    bool perform_extend_out_right(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector);
    void extend_out_left(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, std::size_t max_remaining_non_solid, bool& extension_success);
    void extend_out_right(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, std::size_t max_remaining_non_solid, bool& extension_success);
    bool perform_extend_right_inner_region(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_last_mod_kmer);
    void first_kmer_exhaustive_search(std::string& kmer, const std::vector<std::size_t>& candidates_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t index);
    void first_kmer_consecutive_search(const std::string& kmer, std::vector<C_candidate_path>& candidate_path_vector);
    bool perform_extend_out_first_kmer(const C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector);
    void extend_first_kmer_to_right(C_candidate_path& candidate_path_in, std::vector<C_candidate_path>& candidate_path_vector_out);
    bool check_long_kmers(const std::string& sequence);

    double get_symbol_probability(std::size_t pos);

    void apply_path_to_temporary_read(const C_candidate_path& candidate_path_in, std::string& sequence_out);
    void apply_path_to_read(const C_candidate_path& candidate_path_in, std::string& sequence_out, std::string& qualities_out, std::size_t& num_corrected_errors);
    void apply_path_to_read(const C_candidate_path& candidate_path_in, std::string& sequence_out, std::string& qualities_out, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2);

    void apply_all_paths_to_read(std::vector<std::vector<C_candidate_path> >& candidate_path_vectors, const std::size_t num_non_solid);
    void check_and_extend_first_kmer_to_right(std::vector<C_candidate_path>& candidate_path_vector_first_kmer);

    void check_is_new_potential_indel_5_prime(std::size_t nesting, std::size_t current_symbol_pos, bool last_kmer = false);
    void check_is_new_potential_indel(std::size_t nesting, std::size_t current_symbol_pos, bool last_kmer = false);
    void correct_indel(const std::string& kmer, std::size_t index_kmer, std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid);
    bool correct_last_deletion(const std::string& kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_kmer, std::size_t nesting, std::size_t& checked_changes, bool next_symbol = false);
    void correct_indel_5_prime(const std::string& kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid);
    void correct_indel_3_prime(const std::string& kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid);

    std::vector<C_candidate_path>::iterator choose_best_correction(std::vector<C_candidate_path>& candidate_path_vector);
    void modify_errors(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors);
    void modify_errors_no_rate(C_candidate_path& candidate_path, std::size_t& num_corrected_errors);
    void modify_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2);
    void modify_errors_first_kmer_no_rate(C_candidate_path& candidate_path, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2);
    bool create_modification_path_towards_5_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting);
    bool create_modification_path_with_indels_towards_5_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting);
    bool create_modification_path_towards_3_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t current_kmer_pos, std::size_t correction_nesting);
    bool create_modification_path_with_indels_towards_3_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting);
    bool create_modification_path_towards_3_prime_internal(std::vector<C_candidate_path>& candidate_path_vector, std::size_t current_kmer_pos, std::size_t correction_nesting);
    bool create_modification_path_with_indels_towards_3_prime_internal(std::vector<C_candidate_path>& candidate_path_vector, std::size_t current_kmer_pos, std::size_t correction_nesting);
    bool create_modification_path_towards_3_prime_single_kmer(std::vector<C_candidate_path>& candidate_path_vector, const std::vector<std::size_t>& candidates_indexes, float kmer_quality);
    bool create_modification_path_with_single_substitution(std::vector<C_candidate_path>& candidate_path_vector, float kmer_quality, std::size_t pos, char modification);
    bool create_modification_path_with_single_indel(std::vector<C_candidate_path>& candidate_path_vector, float kmer_quality, std::size_t pos, bool insertion = true, char modification = '0');
    void create_modification_path_empty_no_extend(std::vector<C_candidate_path>& candidate_path_vector, float kmer_quality);
};



//----------------------------------------------------------------------
// Deletion of the unnecessary specializations.
//----------------------------------------------------------------------

template<>
void C_correct_read<false>::check_is_new_potential_indel_5_prime(std::size_t nesting, std::size_t current_symbol_pos, bool last_kmer) = delete;

template<>
void C_correct_read<false>::check_is_new_potential_indel(std::size_t nesting, std::size_t current_symbol_pos, bool last_kmer) = delete;

template<>
void C_correct_read<false>::correct_indel(const std::string& kmer, std::size_t index_kmer, std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid) = delete;

template<>
void C_correct_read<false>::correct_indel_5_prime(const std::string& kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid) = delete;

template<>
void C_correct_read<false>::correct_indel_3_prime(const std::string& kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid) = delete;

template<>
bool C_correct_read<false>::correct_last_deletion(const std::string& kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_kmer, std::size_t nesting, std::size_t& checked_changes, bool next_symbol) = delete;

template<>
bool C_correct_read<false>::create_modification_path_with_indels_towards_5_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting) = delete;

template<>
bool C_correct_read<false>::create_modification_path_with_indels_towards_3_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting) = delete;


//----------------------------------------------------------------------
// Forward declarations for compile error workaround.
//----------------------------------------------------------------------

template<>
void C_correct_read<true>::correct_errors_5_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out);

template<>
void C_correct_read<false>::correct_errors_5_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out);

template<>
void C_correct_read<true>::correct_errors_3_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out);

template<>
void C_correct_read<false>::correct_errors_3_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out);

template<>
void C_correct_read<true>::extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);

template<>
void C_correct_read<false>::extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);

template<>
void C_correct_read<true>::extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);

template<>
void C_correct_read<false>::extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);

template<>
void C_correct_read<false>::extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);

template<>
void C_correct_read<true>::extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid);

template<>
bool C_correct_read<false>::perform_extend_right_inner_region(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_last_mod_kmer);

template<>
bool C_correct_read<true>::perform_extend_right_inner_region(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_last_mod_kmer);

template<>
void C_correct_read<false>::apply_path_to_temporary_read(const C_candidate_path& candidate_path_in, std::string& sequence_out);

template<>
void C_correct_read<true>::apply_path_to_temporary_read(const C_candidate_path& candidate_path_in, std::string& sequence_out);

template<>
void C_correct_read<true>::correct_indel_5_prime(const std::string& kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid);

template<>
void C_correct_read<true>::correct_indel_3_prime(const std::string& kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid);

template<>
bool C_correct_read<true>::perform_extend_out_right(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector);

//----------------------------------------------------------------------
// Adds substitution correction to the path.
//----------------------------------------------------------------------

inline void C_candidate_path::add_substitution(std::size_t pos, char symbol, double symbol_probability) {
    modifications.emplace_back(pos, symbol);
    nucleotides_probability *= symbol_probability;
}



//----------------------------------------------------------------------
// Adds insertion correction to the path.
//----------------------------------------------------------------------

inline void C_candidate_path::add_insertion(std::size_t pos) {
    modifications.emplace_back(pos);
    nucleotides_probability *= INSERTION_RATE_PROB;
    --read_length_change;
    ++num_insertions;
}



//----------------------------------------------------------------------
// Adds deletion correction to the path.
//----------------------------------------------------------------------

inline void C_candidate_path::add_deletion(std::size_t pos, char symbol) {
    modifications.emplace_back(pos, symbol, true);
    nucleotides_probability *= DELETION_RATE_PROB;
    ++read_length_change;
}



//----------------------------------------------------------------------
// Determines erroneous regions and calls their correction.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::correct_errors_in_a_read() {
    //--------------------------------------------------
    // STEP 0-0: find solid k-mers in this read
    //--------------------------------------------------
    // variables
    std::size_t num_kmers(read_length - kmer_length + 1);

    std::vector< std::pair<std::size_t, std::size_t> > solid_regions;

    bool is_solid_kmer_prev(false);

    // find solid regions
    std::pair<std::size_t, std::size_t> new_solid_region;
    for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
        const char* current_kmer = sequence_modified.c_str() + it_kmer;

        // k-mer is solid
        float kmer_quality;
        if (query.query_text(current_kmer, kmer_quality) == true) {
            // start point of a solid region
            if (is_solid_kmer_prev == false) {
                new_solid_region.first = it_kmer;
                is_solid_kmer_prev = true;
            }
        }
        else {
            // end point of a solid region
            if (is_solid_kmer_prev == true) {
                new_solid_region.second = it_kmer - 1;

                if (new_solid_region.second < new_solid_region.first) {
                    c_err << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                solid_regions.push_back(new_solid_region);

                is_solid_kmer_prev = false;
            }
        }
    }


    // last solid region
    if (is_solid_kmer_prev == true) {
        new_solid_region.second = num_kmers - 1;
        solid_regions.push_back(new_solid_region);
    }

    //--------------------------------------------------
    // STEP 0-1: adjust solid regions using quality scores
    //--------------------------------------------------
    // solid_regions_org: indexes that are not modified
    // when the indices are reached, all kinds of modifications (A/C/G/T) should be made and checked
    // at least one solid region
    if (solid_regions.size() > 0) {
        // exceptional case: only one solid island that covers the entire read
        if ((solid_regions.size() != 1) || (solid_regions[0].first != 0) || (solid_regions[0].second != (read_length - kmer_length))) {
            // at least two solid k-mer islands
            if (solid_regions.size() > 1) {
                // check the distance between every two solid k-mer islands
                bool flag_short_distance(false);

                for (std::size_t it_sr = 0; it_sr < solid_regions.size() - 1; it_sr++) {
                    if ((solid_regions[it_sr + 1].first - solid_regions[it_sr].second) < kmer_length) {
                        flag_short_distance = true;
                    }
                }

                if (flag_short_distance == true) {
                    std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;

                    // each solid island
                    for (std::size_t it_sr = 0; it_sr < solid_regions.size(); it_sr++) {
                        // each base in the solid island (0-base)
                        std::size_t num_low_quality_base(0);
                        // set an initial value to avoid a compilation warning
                        std::size_t index_prev_low_quality_base(0);
                        for (size_t it_base = solid_regions[it_sr].first; it_base < solid_regions[it_sr].second + kmer_length; it_base++) {
                            // a current base has a low quality score
                            if (is_low_LUT[quality_score[it_base]]) {
                                num_low_quality_base++;

                                // first low quality base
                                if (num_low_quality_base == 1) {
                                    // the low quality base is not in the first k-mer of the solid island
                                    if (it_base >= (solid_regions[it_sr].first + kmer_length)) {
                                        // add the left most high quality region to a temporary vector
                                        solid_regions_tmp.emplace_back(solid_regions[it_sr].first, (it_base - kmer_length));
                                    }
                                }
                                // not first low quality base
                                else {
                                    if ((it_base - index_prev_low_quality_base) > kmer_length) {
                                        solid_regions_tmp.emplace_back((index_prev_low_quality_base + 1), (it_base - kmer_length));
                                    }
                                }

                                index_prev_low_quality_base = it_base;
                            }
                        }

                        // process the bases to the right of the rightmost low quality base
                        if (num_low_quality_base > 0) {
                            if (solid_regions[it_sr].second >= (index_prev_low_quality_base + kmer_length)) {
                                solid_regions_tmp.emplace_back((index_prev_low_quality_base + kmer_length), solid_regions[it_sr].second);
                            }
                        }
                        // no low quality base
                        // add the current solid island
                        else {
                            solid_regions_tmp.push_back(solid_regions[it_sr]);
                        }
                    }

                    solid_regions.swap(solid_regions_tmp);
                }
            }
            // only one solid k-mer island
            else if (solid_regions.size() == 1) {
                std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;

                std::size_t num_low_quality_base(0);
                std::size_t prev_low_quality_index(0);

                // each base in the solid island (0-base)
                for (size_t it_base = solid_regions[0].first; it_base < solid_regions[0].second + kmer_length; it_base++) {
                    // a current base has a low quality score
                    if (is_low_LUT[quality_score[it_base]]) {
                        num_low_quality_base++;

                        // first low quality base
                        if (num_low_quality_base == 1) {
                            if ((it_base - solid_regions[0].first) >= (kmer_length + MIN_SOLID_LENGTH - 1)) {
                                solid_regions_tmp.emplace_back(solid_regions[0].first, (it_base - kmer_length));
                            }

                            prev_low_quality_index = it_base;
                        }
                        // not first low quality base
                        else {
                            if ((it_base - prev_low_quality_index) >= (kmer_length + MIN_SOLID_LENGTH)) {
                                solid_regions_tmp.emplace_back((prev_low_quality_index + 1), (it_base - kmer_length));
                            }

                            prev_low_quality_index = it_base;
                        }
                    }
                }

                // the above is done only when this procedure does not remove the only solid island
                if (solid_regions_tmp.size() > 0) {
                    solid_regions.swap(solid_regions_tmp);
                }
            }
        }
    }

    //--------------------------------------------------
    // STEP 0-2: remove short solid regions
    //--------------------------------------------------
    if (solid_regions.size() > 0) {
        std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;

        for (std::size_t it_region = 0; it_region < solid_regions.size(); it_region++) {
            if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) >= MIN_SOLID_LENGTH) {
                solid_regions_tmp.push_back(solid_regions[it_region]);
            }
        }

        solid_regions.swap(solid_regions_tmp);
    }

    //--------------------------------------------------
    // STEP 0-3: remove short non-solid regions
    //--------------------------------------------------
    if (solid_regions.size() > 0) {
        std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
        solid_regions_tmp.push_back(solid_regions[0]);

        if (solid_regions.size() > 1) {
            for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
                if ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < MIN_NON_SOLID_LENGTH) {
                    solid_regions_tmp[solid_regions_tmp.size() - 1].second = solid_regions[it_region].second;
                }
                else {
                    solid_regions_tmp.push_back(solid_regions[it_region]);
                }
            }
        }
        solid_regions.swap(solid_regions_tmp);
    }

    //--------------------------------------------------
    // STEP 0-4: reduce the size of solid regions
    //--------------------------------------------------
    //if (solid_regions.size() > 1) {
    //    for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
    //        // (length of a non-solid region < kmer_length) && (length of a non-solid region >= kmer_length - FP_SUSPECT_LENGTH(default: 1))
    //        if (((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < kmer_length) &&
    //                ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) >= kmer_length - FP_SUSPECT_LENGTH)) {
    //            // length of the right solid region > FP_SUSPECT_LENGTH(default: 1)
    //            if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) > FP_SUSPECT_LENGTH) {
    //                solid_regions[it_region].first += FP_SUSPECT_LENGTH;
    //            }
    //
    //            // length of the left solid region > FP_SUSPECT_LENGTH(default: 1)
    //            if ((solid_regions[it_region - 1].second - solid_regions[it_region - 1].first + 1) > FP_SUSPECT_LENGTH) {
    //                solid_regions[it_region - 1].second -= FP_SUSPECT_LENGTH;
    //            }
    //        }
    //    }
    //}

    //--------------------------------------------------
    // STEP 0-5: remove a solid region that makes a non-solid reiong shorter than k
    //--------------------------------------------------
    if (solid_regions.size() == 2) {
        // the first solid region starts from the first k-mer
        if (solid_regions[0].first == 0) {
            // the distance between two regions is shorter than k
            if ((solid_regions[1].first - solid_regions[0].second) < (kmer_length + 1)) {
                // remove the second solid region
                solid_regions.erase(solid_regions.begin() + 1);
            }
        }
            // the second solid region ends in the last k-mer
        else if (solid_regions[1].second == (read_length - kmer_length)) {
            // the distance between two regions is shorter than k
            if ((solid_regions[1].first - solid_regions[0].second) < (kmer_length + 1)) {
                // the length of the second solid region is >= 10% of the sequence length
                if ((solid_regions[1].second - solid_regions[1].first + 1) >= (read_length * 0.1)) {
                    // the length of the first solid region is < 10% of the sequence length
                    if ((solid_regions[0].second - solid_regions[0].first + 1) < (read_length * 0.1)) {
                        // remove the second solid region
                        solid_regions.erase(solid_regions.begin());
                    }
                }
            }
        }
    }

    //--------------------------------------------------
    // STEP 0-6: check the quality scores of right side of each solid k-mer region
    //--------------------------------------------------
    // at least one solid region
    if (solid_regions.size() > 0) {
        // 1 - (n - 1) solid region
        for (std::size_t it_sr = 0; it_sr < (solid_regions.size() - 1); it_sr++) {
            const std::size_t max_adjust = std::min((std::size_t)SOLID_REGION_ADJUST_RANGE, solid_regions[it_sr].second - solid_regions[it_sr].first);
            // sufficient solid regions length
            if ((solid_regions[it_sr].second - solid_regions[it_sr].first) >= max_adjust) {
                for (std::size_t it_adjust = solid_regions[it_sr].second; it_adjust > (solid_regions[it_sr].second - max_adjust); it_adjust--) {
                    // low quality score
                    if (is_low_LUT[quality_score[it_adjust + kmer_length - 1]] ||
                            is_low_LUT[quality_score[it_adjust]]
                            ) {
                        solid_regions[it_sr].second = it_adjust - 1;
                        break;
                    }
                }
            }
        }

        for (std::size_t it_sr = 1; it_sr < solid_regions.size(); it_sr++) {
            const std::size_t max_adjust = std::min((std::size_t)SOLID_REGION_ADJUST_RANGE, solid_regions[it_sr].second - solid_regions[it_sr].first);
            if ((solid_regions[it_sr].second - solid_regions[it_sr].first) >= max_adjust) {
                const std::size_t region_begin = solid_regions[it_sr].first;
                for (std::size_t it_adjust = solid_regions[it_sr].first; it_adjust < (region_begin + max_adjust); it_adjust++) {
                    if (is_low_LUT[quality_score[it_adjust]] ||
                            is_low_LUT[quality_score[it_adjust + kmer_length - 1]]
                            ) {
                        solid_regions[it_sr].first = it_adjust + 1;
                    }
                }
            }
        }

        // last solid region
        std::size_t index_solid_region(solid_regions.size() - 1);

        // non-solid k-mers exist at the 3-prime end
        std::size_t max_adjust = std::min((std::size_t)SOLID_REGION_ADJUST_RANGE, solid_regions[index_solid_region].second - solid_regions[index_solid_region].first);
        // sufficient solid regions length
        if ((solid_regions[index_solid_region].second - solid_regions[index_solid_region].first) >= max_adjust) {
            const std::size_t region_end = solid_regions[index_solid_region].second;
            for (std::size_t it_adjust = solid_regions[index_solid_region].second; it_adjust > (region_end - max_adjust); it_adjust--) {
                // low quality score
                if (((std::size_t)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < QS_CUTOFF) {
                    solid_regions[index_solid_region].second = it_adjust - 1;
                }
            }
        }

        // non-solid k-mers exist at the 5-prime end
        max_adjust = std::min((std::size_t)SOLID_REGION_ADJUST_RANGE, solid_regions[0].second - solid_regions[0].first);
        // sufficient solid regions length
        if ((solid_regions[0].second - solid_regions[0].first) >= max_adjust) {
            const std::size_t region_start = solid_regions[0].first;
            for (std::size_t it_adjust = solid_regions[0].first; it_adjust < (region_start + max_adjust); it_adjust++) {
                // low quality score
                if (is_low_LUT[quality_score[it_adjust]]) {
                    solid_regions[0].first = it_adjust + 1;
                }
            }
        }
    }

    //--------------------------------------------------
    // STEP 0-7: check whether a short non-solid region still exists
    //--------------------------------------------------
    bool short_non_solid_region(false);
    std::size_t short_non_solid_region_length = CORRECT_INDEL ? kmer_length - 1 : kmer_length; // regions of length k-1 are accepted, as they can contain single deletion
    if (solid_regions.size() > 1) {
        for (std::size_t it_sr = 1; it_sr < (solid_regions.size() - 1); it_sr++) {
            if ((solid_regions[it_sr].first - solid_regions[it_sr - 1].second) <= short_non_solid_region_length) {
                short_non_solid_region = true;
                break;
            }
        }
    }

    //--------------------------------------------------
    // correct errors
    //--------------------------------------------------

    if ((solid_regions.size() > 0) && (short_non_solid_region == false)) {
        const bool front_region_exists = solid_regions[0].first > 0;
        const bool back_region_exists = solid_regions[solid_regions.size() - 1].second < (read_length - kmer_length);

        // storage for all potential non solid regions, including back and front ones
        const std::size_t num_non_solid = solid_regions.size() + 1;
        std::vector<std::vector<C_candidate_path> > candidate_path_vectors(num_non_solid);

        //--------------------------------------------------
        // STEP 1-1: Correct errors between solid regions
        //--------------------------------------------------
        if (solid_regions.size() > 1) {
            if (CORRECT_INDEL) {
                // for each solid region
                for (std::size_t it_region = 1; it_region < num_non_solid - 1; it_region++) {
                    if ((((solid_regions[it_region].first - 1) - (solid_regions[it_region - 1].second + 1)) + 1) == kmer_length - 1) {
                        ++num_single_deletions;
                    }
                    if ((((solid_regions[it_region].first - 1) - (solid_regions[it_region - 1].second + 1)) + 1) >= kmer_length - 1) {
                        correct_errors_between_solid_regions(
                            (solid_regions[it_region - 1].second + 1),
                            (solid_regions[it_region].first - 1),
                            candidate_path_vectors[it_region]
                            );
                    }
                }
            }
            else {
                // for each solid region
                for (std::size_t it_region = 1; it_region < num_non_solid - 1; it_region++) {
                    if ((((solid_regions[it_region].first - 1) - (solid_regions[it_region - 1].second + 1)) + 1) >= kmer_length) {
                        correct_errors_between_solid_regions(
                            (solid_regions[it_region - 1].second + 1),
                            (solid_regions[it_region].first - 1),
                            candidate_path_vectors[it_region]
                            );
                    }
                }
            }
        }

        //--------------------------------------------------
        // STEP 1-2: Correct errors in the 5' end
        //--------------------------------------------------
        // the first solid region does not start from the 0-th k-mer in a read
        if (front_region_exists) {
            correct_errors_5_prime_end(solid_regions[0].first - 1, candidate_path_vectors.front());
        }

        //--------------------------------------------------
        // STEP 1-3: Correct errors in the 3' end
        //--------------------------------------------------
        // the last solid region does not end in the last k-mer in a read
        if (back_region_exists) {
            correct_errors_3_prime_end(solid_regions[solid_regions.size() - 1].second + 1, candidate_path_vectors.back());
        }

        // some erroneous regions exist
        if (solid_regions.size() > 1 || front_region_exists || back_region_exists) {
            // check long k-mers (if enabled) and apply changes to the read
            apply_all_paths_to_read(candidate_path_vectors, num_non_solid);
        }
    }
    //--------------------------------------------------
    // no solid region or short weak regions
    //--------------------------------------------------
    else {
        //--------------------------------------------------
        // STEP 2-1: Correct errors in the first k-mer
        //--------------------------------------------------
        // find potentially wrong bases
        std::vector<C_candidate_path> candidate_path_vector_first_kmer;

        correct_errors_first_kmer(candidate_path_vector_first_kmer);

        //--------------------------------------------------
        // STEP 2-2: extend candidate paths to the right
        //--------------------------------------------------
        //--------------------------------------------------
        // STEP 2-3: choose a final one in if possible and optionally verify with long k-mers
        //--------------------------------------------------

        if (candidate_path_vector_first_kmer.size() > 0) {
            check_and_extend_first_kmer_to_right(candidate_path_vector_first_kmer);
        }
    }
}



//----------------------------------------------------------------------
// Init set of variables values for single read correction
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::prepare_corrector() {
    read_length = sequence_modified.length();

    if (read_length > modifications_sequence_stack.size()) {
        modifications_sequence_stack.resize(read_length);
    }
}



//----------------------------------------------------------------------
// Corrects errors in a region situated between correct regions.
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::correct_errors_between_solid_regions(const std::size_t index_start, const std::size_t index_end, std::vector<C_candidate_path>& candidate_path_vector_out) {
    //--------------------------------------------------
    // from i-th region to (i + 1)-th region
    //--------------------------------------------------
    // i-th solid region | non-solid | (i+1)-th solid region
    // --------------------------------------- read
    //              |-----|                    (index_start)-th k-mer
    //                              |-----|    (index_end)-th k-mer: last non-solid k-mer
    //                         |-----|         (index_last_mod)-th k-mer = (index_end - kmer_length + 1)-th k-mer: last k-mer that can be modified
    //                               |----|    This regions should not be modified
    //--------------------------------------------------

    // index of the k-mer that can be modified
    // k-mers that are overlapped with a solid regioin cannot be modified
    std::size_t index_last_mod(index_end - kmer_length + 1);

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    current_region_index_start = (int)index_start;

    if (index_start > index_last_mod) { // can happen, if region contains just a deletion
        max_current_indels = 1;

        std::size_t checked_changes = CHECK_MAX_CHANGES;

        correct_last_deletion(kmer_initial, candidate_path_vector_out, index_start, 0, checked_changes, true);
    }
    else {
        max_current_indels = std::min(MAX_INDELS, static_cast<int>(std::ceil((index_last_mod - index_start + 1) * MAX_INDELS_RATE)));

        const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (index_last_mod - index_start + 1))), MAX_CHANGES_IN_REGION_MIN);

        const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((index_last_mod - index_start + 1) * MAX_NON_SOLID));

        std::size_t checked_changes = CHECK_MAX_CHANGES;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                return;
            }
            --checked_changes;
            // make a change
            kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            float kmer_quality;
            if (query.query_text(kmer_initial, kmer_quality) == true) {
                modifications_sequence_stack[0].modification = NEOCLEOTIDE[it_alter];
                modifications_sequence_stack[0].quality = kmer_quality;

                // if this k-mer is the last k-mer that can be modified
                // running extend_a_kmer_right is not needed any more
                if (index_start == index_last_mod) {
                    // generate a new path
                    if (!create_modification_path_towards_3_prime_internal(candidate_path_vector_out, index_start, 0) || candidate_path_vector_out.empty()) {
                        if (!correct_last_deletion(kmer_initial, candidate_path_vector_out, index_start, 1, checked_changes)) {
                            if (sequence_modified[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
                                check_is_new_potential_indel(0, index_start, true);
                            }
                        }
                    }
                }
                else { // index_start < index_last_mod
                       // trace this kmer recursively and update candidate_path_vector_tmp
                    extend_a_kmer(
                        kmer_initial,
                        index_start + 1,
                        index_last_mod,
                        candidate_path_vector_out,
                        1,
                        checked_changes,
                        NEOCLEOTIDE[it_alter] == sequence_modified[index_start] ? max_remaining_changes : max_remaining_changes - 1,
                        max_remaining_non_solid
                        );
                }
            }

            if (sequence_modified[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
                correct_indel(kmer_initial, index_start, index_last_mod, candidate_path_vector_out, 0, checked_changes, max_remaining_changes, max_remaining_non_solid);
            }
        }

    }
}



template<>
inline void C_correct_read<false>::correct_errors_between_solid_regions(const std::size_t index_start, const std::size_t index_end, std::vector<C_candidate_path>& candidate_path_vector_out) {
    //--------------------------------------------------
    // from i-th region to (i + 1)-th region
    //--------------------------------------------------
    // i-th solid region | non-solid | (i+1)-th solid region
    // --------------------------------------- read
    //              |-----|                    (index_start)-th k-mer
    //                              |-----|    (index_end)-th k-mer: last non-solid k-mer
    //                         |-----|         (index_last_mod)-th k-mer = (index_end - kmer_length + 1)-th k-mer: last k-mer that can be modified
    //                               |----|    This regions should not be modified
    //--------------------------------------------------

    // index of the k-mer that can be modified
    // k-mers that are overlapped with a solid regioin cannot be modified
    std::size_t index_last_mod(index_end - kmer_length + 1);

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (index_last_mod - index_start + 1))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((index_last_mod - index_start + 1) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (checked_changes == 0) {
            return;
        }
        --checked_changes;

        // make a change
        kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query.query_text(kmer_initial, kmer_quality) == true) {
            modifications_sequence_stack[0].modification = NEOCLEOTIDE[it_alter];
            modifications_sequence_stack[0].quality = kmer_quality;

            // if this k-mer is the last k-mer that can be modified
            // running extend_a_kmer_right is not needed any more
            if (index_start == index_last_mod) {
                // generate a new path
                create_modification_path_towards_3_prime_internal(candidate_path_vector_out, index_start, 0);
            }
            else {
                // trace this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer(
                    kmer_initial,
                    index_start + 1,
                    index_last_mod,
                    candidate_path_vector_out,
                    1,
                    checked_changes,
                    NEOCLEOTIDE[it_alter] == sequence_modified[index_start] ? max_remaining_changes : max_remaining_changes - 1,
                    max_remaining_non_solid
                    );
            }
        }
    }
}



//----------------------------------------------------------------------
// Corrects errors situated on a 5' end of the read
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::correct_errors_5_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out) {
    // |  non-solid  | 1st solid region
    // |--------------------------------------| read
    //         |-----|                          (index_start)-th k-mer
    //--------------------------------------------------

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    current_region_index_start = (int)index_start;

    max_current_indels = std::min(MAX_INDELS, static_cast<int>(std::ceil((index_start + 1) * MAX_INDELS_RATE)));

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (index_start + 1))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((index_start + 1) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (checked_changes == 0) {
            return;
        }
        --checked_changes;

        // make a change
        kmer_initial[0] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query.query_text(kmer_initial, kmer_quality) == true) {
            modifications_sequence_stack[0].modification = NEOCLEOTIDE[it_alter];
            modifications_sequence_stack[0].quality = kmer_quality;

            // if this k-mer is the first k-mer in a read
            // running extend_a_kmer_5_prime_end is not needed any more
            if (index_start == 0) {
                create_modification_path_towards_5_prime(candidate_path_vector_out, 0);
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer_5_prime_end(
                    kmer_initial,
                    index_start - 1,
                    candidate_path_vector_out,
                    1,
                    checked_changes,
                    NEOCLEOTIDE[it_alter] == sequence_modified[index_start] ? max_remaining_changes : max_remaining_changes - 1,
                    max_remaining_non_solid
                    );
            }
        }
    }
    correct_indel_5_prime(kmer_initial, index_start, candidate_path_vector_out, 0, checked_changes, max_remaining_changes, max_remaining_non_solid);
}



template<>
inline void C_correct_read<false>::correct_errors_5_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out) {
    // |  non-solid  | 1st solid region
    // |--------------------------------------| read
    //         |-----|                          (index_start)-th k-mer
    //--------------------------------------------------

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (index_start + 1))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((index_start + 1) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (checked_changes == 0) {
            return;
        }
        --checked_changes;

        // make a change
        kmer_initial[0] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query.query_text(kmer_initial, kmer_quality) == true) {
            modifications_sequence_stack[0].modification = NEOCLEOTIDE[it_alter];
            modifications_sequence_stack[0].quality = kmer_quality;

            // if this k-mer is the first k-mer in a read
            // running extend_a_kmer_5_prime_end is not needed any more
            if (index_start == 0) {
                create_modification_path_towards_5_prime(candidate_path_vector_out, 0);
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer_5_prime_end(
                    kmer_initial,
                    index_start - 1,
                    candidate_path_vector_out,
                    1,
                    checked_changes,
                    NEOCLEOTIDE[it_alter] == sequence_modified[index_start] ? max_remaining_changes : max_remaining_changes - 1,
                    max_remaining_non_solid
                    );
            }
        }
    }
}



//----------------------------------------------------------------------
// Corrects errors situated on a 3' end of the read
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::correct_errors_3_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out) {
    //  last solid region | non-solid region |
    // --------------------------------------| read
    //               |-----|                   (index_start)-th k-mer
    //--------------------------------------------------

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    current_region_index_start = (int)index_start;

    max_current_indels = std::min(MAX_INDELS, static_cast<int>(std::ceil((read_length - index_start - kmer_length + 1) * MAX_INDELS_RATE)));

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (read_length - index_start - kmer_length + 1))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((read_length - index_start - kmer_length + 1) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (checked_changes == 0) {
            return;
        }
        --checked_changes;
        // make a change
        kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query.query_text(kmer_initial, kmer_quality) == true) {
            modifications_sequence_stack[0].modification = NEOCLEOTIDE[it_alter];
            modifications_sequence_stack[0].quality = kmer_quality;

            // if this k-mer is the last k-mer in a read
            // running extend_a_kmer_3_prime_end is not needed any more
            if (index_start == (read_length - kmer_length)) {
                create_modification_path_towards_3_prime(candidate_path_vector_out, index_start, 0);
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_out
                extend_a_kmer_3_prime_end(
                    kmer_initial,
                    index_start + 1,
                    candidate_path_vector_out,
                    1,
                    checked_changes,
                    NEOCLEOTIDE[it_alter] == sequence_modified[index_start + kmer_length - 1] ? max_remaining_changes : max_remaining_changes - 1,
                    max_remaining_non_solid
                    );
            }
        }
    }
    correct_indel_3_prime(kmer_initial, index_start, candidate_path_vector_out, 0, checked_changes, max_remaining_changes, max_remaining_non_solid);
}



template<>
inline void C_correct_read<false>::correct_errors_3_prime_end(const std::size_t index_start, std::vector<C_candidate_path>& candidate_path_vector_out) {
    //  last solid region | non-solid region |
    // --------------------------------------| read
    //               |-----|                   (index_start)-th k-mer
    //--------------------------------------------------

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (read_length - index_start - kmer_length + 1))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((read_length - index_start - kmer_length + 1) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;
    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (checked_changes == 0) {
            return;
        }
        --checked_changes;

        // make a change
        kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query.query_text(kmer_initial, kmer_quality) == true) {
            modifications_sequence_stack[0].modification = NEOCLEOTIDE[it_alter];
            modifications_sequence_stack[0].quality = kmer_quality;

            // if this k-mer is the last k-mer in a read
            // running extend_a_kmer_3_prime_end is not needed any more
            if (index_start == (read_length - kmer_length)) {
                create_modification_path_towards_3_prime(candidate_path_vector_out, index_start, 0);
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_out
                extend_a_kmer_3_prime_end(
                    kmer_initial,
                    index_start + 1,
                    candidate_path_vector_out,
                    1,
                    checked_changes,
                    NEOCLEOTIDE[it_alter] == sequence_modified[index_start + kmer_length - 1] ? max_remaining_changes : max_remaining_changes - 1,
                    max_remaining_non_solid
                    );
            }
        }
    }
}



//----------------------------------------------------------------------
// Checks if the sequence is covered by long k-mers.
//----------------------------------------------------------------------
template<bool CORRECT_INDEL>
bool C_correct_read<CORRECT_INDEL>::check_long_kmers(const std::string& sequence) {
    if (!query.use_long_kmer()) {
        return true;
    }

    const char* kmer = sequence.c_str();
    for (std::size_t it_base = 0; it_base < sequence.length() - long_kmer_length + 1; ++it_base, ++kmer) {
        float tmp;
        if (!query.query_text_long_kmer(kmer, tmp)) {
            return false;
        }
    }
    return true;
}



//----------------------------------------------------------------------
// Corrects errors in the first k-mer of the read.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline void C_correct_read<CORRECT_INDEL>::correct_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector_out) {
    std::string first_kmer(sequence_modified.substr(0, kmer_length));
    ++num_first_kmer_corrections;

    std::vector<std::size_t> low_qs_indexes;

    for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
        if (is_low_LUT[quality_score[it_bases]]) {
            low_qs_indexes.push_back(it_bases);
        }
    }

    // correct errors if the number of low-quality bases is smaller than the threshold
    if ((low_qs_indexes.size() <= MAX_LOW_QS_BASES) && (low_qs_indexes.size() > 0)) {
        std::string kmer(first_kmer);
        first_kmer_exhaustive_search(
                kmer,
                low_qs_indexes,
                candidate_path_vector_out,
                0
                );

        // no candidate path is found
        if (candidate_path_vector_out.size() == 0) {
            first_kmer_consecutive_search(first_kmer, candidate_path_vector_out);
        }
    }
        // no low-quality base or too many low-quality bases
    else {
        float kmer_quality;
        if (query.query_text(first_kmer, kmer_quality) == true) {
            create_modification_path_empty_no_extend(candidate_path_vector_out, kmer_quality);
        }
        else {
            if (low_qs_indexes.size() == 0) {
                first_kmer_consecutive_search(first_kmer, candidate_path_vector_out);
            }
            else {
                // (quality, position)
                std::vector<std::pair<char, std::size_t> > qualities_vector;
                qualities_vector.reserve(kmer_length);

                for (std::size_t it_qualities = 0; it_qualities < kmer_length; it_qualities++) {
                    qualities_vector.emplace_back(quality_score[it_qualities], it_qualities);
                }
                std::sort(qualities_vector.begin(), qualities_vector.end());

                std::vector<std::size_t> candidates;

                for (std::size_t it = 0; it < MAX_LOW_QS_BASES; ++it) {
                    candidates.push_back(qualities_vector[it].second);
                }

                std::string kmer(first_kmer);
                first_kmer_exhaustive_search(
                    kmer,
                    candidates,
                    candidate_path_vector_out,
                    0
                    );
                // no candidate path is found
                if (candidate_path_vector_out.size() == 0) {
                    first_kmer_consecutive_search(first_kmer, candidate_path_vector_out);
                }
            }
        }
    }

    if (CORRECT_INDEL) {
        if (candidate_path_vector_out.empty()) {
            // try insertion
            // check every position except the first one
            for (std::size_t it = 1; it < kmer_length; ++it) {
                std::string kmer_tmp(sequence_modified.substr(0, kmer_length + 1));
                kmer_tmp.erase(it, 1);

                float kmer_quality;
                if (query.query_text(kmer_tmp, kmer_quality) == true) {
                    create_modification_path_with_single_indel(candidate_path_vector_out, kmer_quality, it);
                }
            }

            // try deletion
            // check positions between the all neighbours
            for (std::size_t it = 1; it < kmer_length; ++it) {
                // each alternative neocletide
                for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
                    std::string kmer_tmp(sequence_modified.substr(0, kmer_length - 1));
                    kmer_tmp.insert(it, 1, NEOCLEOTIDE[it_alter]);

                    float kmer_quality;
                    if (query.query_text(kmer_tmp, kmer_quality) == true) {
                        create_modification_path_with_single_indel(candidate_path_vector_out, kmer_quality, it, false, NEOCLEOTIDE[it_alter]);
                    }
                }
            }
        }
    }

    if (candidate_path_vector_out.size() > 0) {
        ++num_first_kmer_successes;
    }

    if (candidate_path_vector_out.size() > MAX_FIRST_KMER_POSSIBILITIES) {
        std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;
        candidate_path_vector_tmp_tmp.swap(candidate_path_vector_out);

        // (rate, index)
        std::vector<std::pair<double, std::vector<C_candidate_path>::iterator> > rates;
        rates.reserve(candidate_path_vector_tmp_tmp.size());

        // each candidate path
        for(std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); ++it_path) {
            // each modification
            rates.emplace_back(it_path->get_path_rate(), it_path);
        }

        std::sort(rates.begin(), rates.end());

        for (std::size_t it = 0; it < (std::size_t)MAX_FIRST_KMER_POSSIBILITIES; ++it) {
            candidate_path_vector_out.push_back(*rates[rates.size() - it - 1].second);
        }
    }
}



//----------------------------------------------------------------------
// Tries to correct k-mer by changing one symbol.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline void C_correct_read<CORRECT_INDEL>::first_kmer_exhaustive_search(std::string& kmer, const std::vector<std::size_t>& candidates_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t index) {
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (candidate_path_vector.size() >= MAX_FIRST_KMER_CORRECTION_PATHS) {
            return;
        }
        // make a new k-mer
        kmer[candidates_indexes[index]] = NEOCLEOTIDE[it_alter];

        modifications_sequence_stack[index].modification = NEOCLEOTIDE[it_alter];

        if (index == candidates_indexes.size() - 1) {
            float kmer_quality;
            if (query.query_text(kmer, kmer_quality) == true) {
                create_modification_path_towards_3_prime_single_kmer(candidate_path_vector, candidates_indexes, kmer_quality);
            }
        }
        else {
            first_kmer_exhaustive_search(
                    kmer,
                    candidates_indexes,
                    candidate_path_vector,
                    index + 1
                    );
        }
    }
}



//----------------------------------------------------------------------
// Tries to correct k-mer by changing one symbol.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline void C_correct_read<CORRECT_INDEL>::first_kmer_consecutive_search(const std::string& kmer, std::vector<C_candidate_path>& candidate_path_vector) {
    std::string kmer_tmp(kmer);
    for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // not equal to the original character
            if (kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                // generate a new k-mer
                kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                // add kmer_tmp to candidate_path_tmp if it is solid
                float kmer_quality;
                if (query.query_text(kmer_tmp, kmer_quality) == true) {
                    // generate a new candidate path
                    create_modification_path_with_single_substitution(candidate_path_vector, kmer_quality, it_bases, NEOCLEOTIDE[it_alter]);
                }
            }
        }
        kmer_tmp[it_bases] = kmer[it_bases];
    }
}




//----------------------------------------------------------------------
// Checks if modified first k-mer can be extended.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::perform_extend_out_first_kmer(const C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector) {
    // no modified path
    if (candidate_path.modifications.size() == 0) {
        candidate_path_vector.push_back(candidate_path);
        return false;
    }

    // check the index of the first modified base
    // extension is needed
    //else if (candidate_path_vector_tmp[it_candidates].modifications[0].first < (MAX_EXTENSION - 1)) {
    if (candidate_path.modifications[0].pos < (kmer_length - 1)) {
        // index_smallest_modified
        std::size_t index_smallest_modified(candidate_path.modifications[0].pos);

        // number of bases that should be extended
        std::size_t extend_amount;

        // applied the modified bases to first_kmer
        std::string first_kmer(sequence_modified.substr(0, kmer_length));
        apply_path_to_temporary_read(candidate_path, first_kmer);

        // determine the number of extensions
        // extension amount = kmer_length - index_smallest_modified - 1
        // kmer_length = 11, max_extension = 5
        // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
        // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
        //           |<------------------->|      k = 11
        //           |--------------------------- read
        //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
        //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
        if (index_smallest_modified >= kmer_length - max_extension - 1) {
            extend_amount = kmer_length - index_smallest_modified - 1;
        }
            // kmer_length = 11, max_extension = 5
            // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
            // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
            //           |<------------------->|      k = 11
            //           |--------------------------- read
            //           |<------->|                  index_smallest_modified < 5
        else {
            extend_amount = max_extension;
        }


        bool extension_success(false);
        std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil(extend_amount * MAX_NON_SOLID));

        // generate an initial k-mer
        std::string kmer_initial(first_kmer.substr(0, kmer_length - 1));
        kmer_initial = '0' + kmer_initial;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[0] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            float kmer_quality = 0.0f;
            bool solid_kmer = query.query_text(kmer_initial, kmer_quality);
            if (solid_kmer || max_remaining_non_solid > 0) {
                // if extend_amount == 1
                // running extend_out_left is not needed any more
                if (extend_amount == 1) {
                    if (solid_kmer) {
                        extension_success = true;
                    }
                    break;
                }
                else if (!extension_success) {
                    // trace  this kmer recursively and update candidate_path_vector_tmp
                    extend_out_left(
                            kmer_initial,
                            1,
                            extend_amount,
                            solid_kmer ? max_remaining_non_solid : max_remaining_non_solid - 1,
                            extension_success
                            );
                }
            }
        }

        if (extension_success == true) {
            candidate_path_vector.push_back(candidate_path);
        }
        return extension_success;
    }
    // extension is not needed
    else {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
}



//----------------------------------------------------------------------
// Proceeds correction from modified first k-mer towards 3' end.
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::extend_first_kmer_to_right(C_candidate_path& candidate_path_in, std::vector<C_candidate_path>& candidate_path_vector_out) {
    // generate the first k-mer

    // first k-mer correction can cause situation, when first "k-mer" has another length
    const int read_length_change = candidate_path_in.get_read_length_change();

    std::string first_kmer(sequence_modified.substr(0, kmer_length));
    int second_kmer_pos = 1;
    if (read_length_change < 0) { // there were an insertion
        if (sequence_modified.length() >= kmer_length + 1) {
            first_kmer = first_kmer.substr(1, kmer_length - 1);
            first_kmer += sequence_modified[kmer_length];
            ++second_kmer_pos;
        }
        // there is no possibility to extend the k-mer
        else {
            perform_extend_out_right(candidate_path_in, candidate_path_vector_out);
            return;
        }
    }

    apply_path_to_temporary_read(candidate_path_in, first_kmer);

    max_current_indels = std::min(MAX_INDELS, static_cast<int>(std::ceil((read_length - kmer_length - second_kmer_pos + 1) * MAX_INDELS_RATE)));

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (read_length - kmer_length - second_kmer_pos + 1))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((read_length - kmer_length - second_kmer_pos + 1) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;

    current_region_index_start = second_kmer_pos;

    extend_a_kmer_3_prime_end(
        first_kmer,
        second_kmer_pos,
        candidate_path_vector_out,
        0,
        checked_changes,
        max_remaining_changes,
        max_remaining_non_solid
        );

    // complete the paths with a modifications of the first k-mer
    for (C_candidate_path& path : candidate_path_vector_out) {
        path.attach_modifications_at_beginning(candidate_path_in);
    }
}



template<>
inline void C_correct_read<false>::extend_first_kmer_to_right(C_candidate_path& candidate_path_in, std::vector<C_candidate_path>& candidate_path_vector_out) {
    // generate the first k-mer
    std::string first_kmer(sequence_modified.substr(0, kmer_length + 1));

    apply_path_to_temporary_read(candidate_path_in, first_kmer);

    current_region_index_start = 1;

    const int max_remaining_changes = std::min(static_cast<int>(std::ceil(MAX_CHANGES_IN_REGION_RATIO * (read_length - kmer_length))), MAX_CHANGES_IN_REGION_MIN);

    const std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil((read_length - kmer_length) * MAX_NON_SOLID));

    std::size_t checked_changes = CHECK_MAX_CHANGES;

    extend_a_kmer_3_prime_end(
        first_kmer,
        1,
        candidate_path_vector_out,
        0,
        checked_changes,
        max_remaining_changes,
        max_remaining_non_solid
        );

    // complete the paths with a modifications of the first k-mer
    for (C_candidate_path& path : candidate_path_vector_out) {
        path.attach_modifications_at_beginning(candidate_path_in);
    }
}



//----------------------------------------------------------------------
// Applies all regions paths to the read.
// Sometimes checks, if some correction path is covered by long k-mers.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::apply_all_paths_to_read(std::vector<std::vector<C_candidate_path> >& candidate_path_vectors, const std::size_t num_non_solid) {
    const bool front_region_exists = !candidate_path_vectors.front().empty();
    const bool back_region_exists = !candidate_path_vectors.back().empty();

    if (query.use_long_kmer()) {
        ++num_checked_with_long_kmer;

        // prepare to check long k-mers
        for (std::size_t it = 0; it < candidate_path_vectors.size(); ++it) {
            std::sort(candidate_path_vectors[it].begin(), candidate_path_vectors[it].end(),
                [](const C_candidate_path& a, const C_candidate_path& b) -> bool {
                return a.get_path_rate() > b.get_path_rate();
            }
            );
        }

        // regions iterators, initially pointing to first elements
        std::vector<std::size_t> regions_paths_iterators(num_non_solid, 0);

        std::string sequence_tmp(sequence_modified);

        double rates_sum = 0.0;
        // apply the best paths
        // traverse with reversed order to workaround changing the read length
        if (back_region_exists) {
            apply_path_to_temporary_read(candidate_path_vectors.back().front(), sequence_tmp);
            rates_sum += candidate_path_vectors.back().front().get_path_rate();
        }
        for (int it = static_cast<int>(num_non_solid) - 2; it >= 1; --it) {
            // can be false if front or back region does not exists, if there are no generated paths, or the region was to short to correct
            if (!candidate_path_vectors[it].empty()) {
                apply_path_to_temporary_read(candidate_path_vectors[it].front(), sequence_tmp);
                rates_sum += candidate_path_vectors[it].front().get_path_rate();
            }
        }
        if (front_region_exists) {
            apply_path_to_temporary_read(candidate_path_vectors.front().front(), sequence_tmp);
            rates_sum += candidate_path_vectors.front().front().get_path_rate();
        }

        bool something_worsen = false;
        // when false - no path can be worsen
        bool found_path_to_worsen = true;
        while (found_path_to_worsen && !check_long_kmers(sequence_tmp)) {
            something_worsen = true;
            found_path_to_worsen = false;
            std::size_t min_diff_region_pos = 0;

            double min_diff = rates_sum; // choose maximal possible value
            for (std::size_t region_it = 0; region_it < num_non_solid; ++region_it) {
                // current region's path
                std::size_t current_path_iterator = regions_paths_iterators[region_it];

                // check if a worse path exists
                // in case of non-existing paths for the region, it will also omit them
                if (candidate_path_vectors[region_it].size() > current_path_iterator + 1) {
                    double diff = candidate_path_vectors[region_it][current_path_iterator].get_path_rate() - candidate_path_vectors[region_it][current_path_iterator + 1].get_path_rate();
                    if (diff < min_diff) {
                        found_path_to_worsen = true;
                        min_diff = diff;
                        min_diff_region_pos = region_it;
                    }
                }
            }
            if (found_path_to_worsen) {
                // choose worsen path
                ++regions_paths_iterators[min_diff_region_pos];

                // apply worsened paths to the read
                sequence_tmp = sequence_modified;
                for (int it = static_cast<int>(num_non_solid) - 1; it >= 0; --it) {
                    if (!candidate_path_vectors[it].empty()) {
                        std::size_t current_path_iterator = regions_paths_iterators[it];
                        apply_path_to_temporary_read(candidate_path_vectors[it][current_path_iterator], sequence_tmp);
                    }
                }
            }
        }

        // read with the applied paths is covered by long k-mers
        if (found_path_to_worsen) {
            if (!something_worsen) {
                ++num_instantly_accepted_with_long_kmer;
            }
            if (back_region_exists) {
                modify_errors_no_rate(candidate_path_vectors.back()[regions_paths_iterators.back()], num_corrected_errors_step1_3);
            }
            for (std::size_t it = num_non_solid - 2; it >= 1; --it) {
                if (!candidate_path_vectors[it].empty()) {
                    modify_errors_no_rate(candidate_path_vectors[it][regions_paths_iterators[it]], num_corrected_errors_step1_1);
                }
            }
            if (front_region_exists) {
                modify_errors_no_rate(candidate_path_vectors.front()[regions_paths_iterators.front()], num_corrected_errors_step1_2);
            }
            if (CORRECT_INDEL) {
                read_length = sequence_modified.length();
            }
        }
        else {
            ++num_filtered_with_long_kmer;
            if (accept_filtered_with_long_kmers) {
                if (back_region_exists && !candidate_path_vectors.back().empty()) {
                    modify_errors_no_rate(candidate_path_vectors.back().front(), num_corrected_errors_step1_3);
                }
                for (std::size_t it = num_non_solid - 2; it >= 1; --it) {
                    if (!candidate_path_vectors[it].empty()) {
                        modify_errors_no_rate(candidate_path_vectors[it].front(), num_corrected_errors_step1_1);
                    }
                }
                if (front_region_exists && !candidate_path_vectors.front().empty()) {
                    modify_errors_no_rate(candidate_path_vectors.front().front(), num_corrected_errors_step1_2);
                }
                if (CORRECT_INDEL) {
                    read_length = sequence_modified.length();
                }
            }
        }
    }
    else {
        // do not use long k-mers
        if (back_region_exists) {
            modify_errors(candidate_path_vectors.back(), num_corrected_errors_step1_3);
        }
        for (std::size_t it = num_non_solid - 2; it >= 1; --it) {
            modify_errors(candidate_path_vectors[it], num_corrected_errors_step1_1);
        }
        if (front_region_exists) {
            modify_errors(candidate_path_vectors.front(), num_corrected_errors_step1_2);
        }
        if (CORRECT_INDEL) {
            read_length = sequence_modified.length();
        }
    }
}



//----------------------------------------------------------------------
// Calls first k-mer extending to right.
// Sometimes checks, if some correction path is covered by long k-mers.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::check_and_extend_first_kmer_to_right(std::vector<C_candidate_path>& candidate_path_vector_first_kmer) {
    if (query.use_long_kmer()) {
        // Actually, the two-dimensional array is not necessary. One dimension would be sufficient, but here we have separation of first and consecutive k-mers corrections.
        // results of the 2-2 step
        std::vector<std::vector<C_candidate_path> > candidate_paths_vector_extending_regions(candidate_path_vector_first_kmer.size());

        // each path
        for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_first_kmer.size(); it_candidates++) {
            extend_first_kmer_to_right(
                candidate_path_vector_first_kmer[it_candidates],
                candidate_paths_vector_extending_regions[it_candidates]
                );
        }

        //--------------------------------------------------
        // STEP 2-3: choose a final one in candidate_paths_vector_extending_regions if possible and verify with long k-mers
        //--------------------------------------------------

        // pair: first k-mer path number, its extending path number
        std::vector<std::pair<std::size_t, std::size_t> > paths_indexes;
        for (std::size_t it_first = 0; it_first < candidate_path_vector_first_kmer.size(); ++it_first) {
            for (std::size_t it_extending = 0; it_extending < candidate_paths_vector_extending_regions[it_first].size(); ++it_extending) {
                paths_indexes.emplace_back(it_first, it_extending);
            }
        }

        if (!paths_indexes.empty()) {
            std::sort(paths_indexes.begin(), paths_indexes.end(), [&candidate_paths_vector_extending_regions](std::pair<std::size_t, std::size_t> a, std::pair<std::size_t, std::size_t> b) -> bool {
                return candidate_paths_vector_extending_regions[a.first][a.second].get_path_rate() > candidate_paths_vector_extending_regions[b.first][b.second].get_path_rate();
            });

            std::size_t it_pairs = 0;
            std::string sequence_tmp;
            bool is_covered_with_long = false;
            do {
                sequence_tmp = sequence_modified;
                const std::size_t first_path_ix = paths_indexes[it_pairs].first;
                const std::size_t extending_path_ix = paths_indexes[it_pairs].second;

                apply_path_to_temporary_read(candidate_paths_vector_extending_regions[first_path_ix][extending_path_ix], sequence_tmp);
                is_covered_with_long = check_long_kmers(sequence_tmp);
            } while (!is_covered_with_long && ++it_pairs < paths_indexes.size());

            if (is_covered_with_long) {
                const std::size_t first_path_ix = paths_indexes[it_pairs].first;
                const std::size_t extending_path_ix = paths_indexes[it_pairs].second;

                modify_errors_first_kmer_no_rate(candidate_paths_vector_extending_regions[first_path_ix][extending_path_ix], num_corrected_errors_step2_1, num_corrected_errors_step2_2);
            }
        }
    }
    else {
        // results of the 2-2 step
        std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

        // each path
        for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_first_kmer.size(); it_candidates++) {
            std::vector<C_candidate_path> candidate_path_vector_first_extend;

            extend_first_kmer_to_right(
                candidate_path_vector_first_kmer[it_candidates],
                candidate_path_vector_first_extend
                );

            //--------------------------------------------------
            // STEP 2-3: choose a final one in candidate_path_vector_tmp_tmp if possible
            //--------------------------------------------------

            std::vector<C_candidate_path>::iterator best_it_path = choose_best_correction(candidate_path_vector_first_extend);

            if (best_it_path != candidate_path_vector_first_extend.end()) {
                candidate_path_vector_tmp_tmp.push_back(*best_it_path);
            }

            // candidiate_path_vector_tmp_tmp: successfully right extended candidates
        }
        modify_errors_first_kmer(candidate_path_vector_tmp_tmp, num_corrected_errors_step2_1, num_corrected_errors_step2_2);
    }
}



//----------------------------------------------------------------------
// Determines if a potential indel occurred within the previous (i.e. located on the right) symbols.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::check_is_new_potential_indel_5_prime(std::size_t nesting, std::size_t current_symbol_pos, bool last_kmer) {
    if (max_current_indels == 0) {
        return;
    }

    // short region can be tested entirely
    if (last_kmer && nesting + 1 < MAX_CHECK_INDEL_NON_SOLID) {
        std::size_t modifications = 0;
        std::size_t indel_pos_result = 0;
        for (int it = static_cast<int>(nesting); it >= 0; --it) {
            if (modifications_sequence_stack[it].insertion || modifications_sequence_stack[it].deletion) {
                // omit rest of the sequence, however do not exclude possibility of indel presence
                break;
            }
            if (modifications_sequence_stack[it].modification != sequence_modified[current_symbol_pos + nesting - it]) {
                ++modifications;
                indel_pos_result = it;
            }
        }
        // in such a short region, the number of modifications suggesting the indel presence can be smaller
        if (modifications >= std::ceil(MAX_CHECK_INDEL_NON_SOLID * (nesting + 1) / (double)MAX_CHECK_INDEL_SYMBOLS)) {
            for (int it = static_cast<int>(nesting); it >= static_cast<int>(indel_pos_result); --it) {
                modifications_sequence_stack[it].potential_indel = true;
            }
            // update number of currently_proposed_indels
            if (currently_proposed_indels > 0) {
                currently_proposed_indels = 0;
                for (int it = static_cast<int>(indel_pos_result) - 1; it >= 0; --it) {
                    if (modifications_sequence_stack[it].potential_indel) {
                        ++currently_proposed_indels;
                    }
                }
            }
            currently_proposed_indels += static_cast<int>(nesting) - static_cast<int>(indel_pos_result) + 1;
        }
        return;
    }

    if (currently_proposed_indels > 0) {
        return;
    }

  
    if (nesting + 1 >= MAX_CHECK_INDEL_NON_SOLID) {
        const std::size_t max_check = std::min(nesting, (std::size_t)MAX_CHECK_INDEL_SYMBOLS);

        std::size_t modifications = 0;
        std::size_t indel_pos_result = 0;

        for (int it = static_cast<int>(nesting); it >= static_cast<int>(nesting - max_check); --it) {
            if (modifications_sequence_stack[it].insertion || modifications_sequence_stack[it].deletion) {
                // omit rest of the sequence, however do not exclude possibility of indel presence
                break;
            }
            if (modifications_sequence_stack[it].modification != sequence_modified[current_symbol_pos + nesting - it]) {
                ++modifications;
                indel_pos_result = it;
            }
        }

        if (modifications >= MAX_CHECK_INDEL_NON_SOLID) {
            assert(!(modifications_sequence_stack[indel_pos_result].insertion || modifications_sequence_stack[indel_pos_result].deletion));
            modifications_sequence_stack[indel_pos_result].potential_indel = true;
            ++currently_proposed_indels;
        }
    }
}



//----------------------------------------------------------------------
// Determines if a potential indel occurred within the previous symbols.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::check_is_new_potential_indel(std::size_t nesting, std::size_t current_symbol_pos, bool last_kmer) {
    if (max_current_indels == 0) {
        return;
    }

    // short region can be tested entirely
    if (last_kmer && nesting + 1 < MAX_CHECK_INDEL_NON_SOLID) {
        std::size_t modifications = 0;
        std::size_t indel_pos_result = 0;
        for (int it = static_cast<int>(nesting); it >= 0; --it) {
            if (modifications_sequence_stack[it].insertion || modifications_sequence_stack[it].deletion) {
                // omit rest of the sequence, however do not exclude possibility of indel presence
                break;
            }
            if (modifications_sequence_stack[it].modification != sequence_modified[current_symbol_pos - nesting + it]) {
                ++modifications;
                indel_pos_result = it;
            }
        }
        // in such a short region, the number of modifications suggesting the indel presence can be smaller
        if (modifications >= std::ceil(MAX_CHECK_INDEL_NON_SOLID * (nesting + 1) / (double)MAX_CHECK_INDEL_SYMBOLS)) {
            for (int it = static_cast<int>(nesting); it >= static_cast<int>(indel_pos_result); --it) {
                modifications_sequence_stack[it].potential_indel = true;
            }
            // update number of currently_proposed_indels
            if (currently_proposed_indels > 0) {
                currently_proposed_indels = 0;
                for (int it = static_cast<int>(indel_pos_result) - 1; it >= 0; --it) {
                    if (modifications_sequence_stack[it].potential_indel) {
                        ++currently_proposed_indels;
                    }
                }
            }
            currently_proposed_indels += static_cast<int>(nesting) - static_cast<int>(indel_pos_result) + 1;
        }
        return;
    }

    if (currently_proposed_indels > 0) {
        return;
    }

    if (nesting + 1 >= MAX_CHECK_INDEL_NON_SOLID) {
        const std::size_t max_check = std::min(nesting, (std::size_t)MAX_CHECK_INDEL_SYMBOLS);

        std::size_t modifications = 0;
        std::size_t indel_pos_result = 0;

        for (int it = static_cast<int>(nesting); it >= static_cast<int>(nesting - max_check); --it) {
            if (modifications_sequence_stack[it].insertion || modifications_sequence_stack[it].deletion) {
                // omit rest of the sequence, however do not exclude possibility of indel presence
                break;
            }
            if (modifications_sequence_stack[it].modification != sequence_modified[current_symbol_pos - nesting + it]) {
                ++modifications;
                indel_pos_result = it;
            }
        }

        if (modifications >= MAX_CHECK_INDEL_NON_SOLID) {
            assert(!(modifications_sequence_stack[indel_pos_result].insertion || modifications_sequence_stack[indel_pos_result].deletion));
            modifications_sequence_stack[indel_pos_result].potential_indel = true;
            ++currently_proposed_indels;
        }
    }
}



//----------------------------------------------------------------------
// Corrects indel and traces recursively.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::correct_indel(const std::string& kmer, std::size_t index_kmer, std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid) {
    if (!modifications_sequence_stack[nesting].potential_indel) {
        return;
    }

    modifications_sequence_stack[nesting].potential_indel = false;
    --currently_proposed_indels;

    if (max_current_indels == 0) {
        return;
    }

    assert(!(modifications_sequence_stack[nesting].insertion || modifications_sequence_stack[nesting].deletion));
    if (nesting > 0 && (modifications_sequence_stack[nesting - 1].insertion || modifications_sequence_stack[nesting - 1].deletion)) {
        return;
    }

    --max_current_indels;

    bool was_insertion = false, was_deletion = false;
    if (nesting >= MIN_INS_DEL_DIST) {
        for (int i = static_cast<int>(nesting) - 1; i >= static_cast<int>(nesting) - MIN_INS_DEL_DIST; --i) {
            if (modifications_sequence_stack[i].insertion) {
                was_insertion = true;
            }
            else if (modifications_sequence_stack[i].deletion) {
                was_deletion = true;
            }
        }
    }


    // try insertion

    if (!was_deletion) {
        // revert potential substitution correction
        modifications_sequence_stack[nesting].modification = sequence_modified[index_kmer + kmer_length - 1];
        modifications_sequence_stack[nesting].insertion = true;

        if (index_kmer == index_last_mod) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                modifications_sequence_stack[nesting].insertion = false;
                return;
            }
            create_modification_path_with_indels_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting);
        }
        else {
            if (checked_changes == 0) {
                modifications_sequence_stack[nesting].insertion = false;
                return;
            }
            --checked_changes;
            std::string insertion_kmer("0"); // dummy symbol, will be removed by extend_a_kmer function
            insertion_kmer += kmer.substr(0, kmer_length - 1);

            extend_a_kmer(
                insertion_kmer,
                index_kmer + 1,
                index_last_mod,
                candidate_path_vector,
                nesting + 1,
                checked_changes,
                max_remaining_changes,
                max_remaining_non_solid
                );
        }

        modifications_sequence_stack[nesting].insertion = false;
    }

    // try deletion

    if (!was_insertion) {
        std::string deletion_kmer(kmer);
        modifications_sequence_stack[nesting].deletion = true;
        // try different symbols for deletion
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                modifications_sequence_stack[nesting].deletion = false;
                return;
            }

            deletion_kmer.back() = NEOCLEOTIDE[it_alter];

            float kmer_quality;
            if (query.query_text(deletion_kmer, kmer_quality)) {
                modifications_sequence_stack[nesting].quality = kmer_quality;
                modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                --checked_changes;
                extend_a_kmer(
                    deletion_kmer,
                    index_kmer, // next k-mer and get again the current one
                    index_last_mod,
                    candidate_path_vector,
                    nesting + 1,
                    checked_changes,
                    max_remaining_changes,
                    max_remaining_non_solid
                    );
            }
        }
        modifications_sequence_stack[nesting].deletion = false;
    }

    ++max_current_indels;
}



//----------------------------------------------------------------------
// Corrects indel and does not trace recursively.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::correct_last_deletion(const std::string& kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_kmer, std::size_t nesting, std::size_t& checked_changes, bool prev_symbol) {
    // If the following conditions break function, it returns true, as it means, that
    // continuing of indels in the current neighbourhood is aimless.
    if (nesting > 0 && (modifications_sequence_stack[nesting - 1].insertion || modifications_sequence_stack[nesting - 1].deletion)) {
        return true;
    }

    assert(max_current_indels >= 0);
    if (max_current_indels == 0) {
        return true;
    }

    bool extension_success = false;
    --max_current_indels;

    std::string deletion_kmer(kmer);
    if (!prev_symbol) {
        // Move one position ahead, as we correct a deletion after the last symbol in the region
        deletion_kmer = deletion_kmer.substr(1);
        deletion_kmer += '0';
    }
    modifications_sequence_stack[nesting].deletion = true;
    // try different symbols for deletion
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        if (checked_changes == 0) {
            modifications_sequence_stack[nesting].deletion = false;
            return true; // signals success, but returning true would prevent checking indels presence
        }
        --checked_changes;

        deletion_kmer.back() = NEOCLEOTIDE[it_alter];

        float kmer_quality;
        if (query.query_text(deletion_kmer, kmer_quality)) {
            modifications_sequence_stack[nesting].quality = kmer_quality;
            modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                modifications_sequence_stack[nesting].deletion = false;
                return true;
            }
            if (create_modification_path_with_indels_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting)) {
                extension_success = true;
            }
        }
    }
    modifications_sequence_stack[nesting].deletion = false;
    ++max_current_indels;

    return extension_success;
}

//----------------------------------------------------------------------
// Corrects indel situated on a 5' end of the read and traces recursively.
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::correct_indel_5_prime(const std::string & kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid) {
    if (!modifications_sequence_stack[nesting].potential_indel) {
        return;
    }

    modifications_sequence_stack[nesting].potential_indel = false;
    --currently_proposed_indels;

    if (max_current_indels == 0) {
        return;
    }

    assert(!(modifications_sequence_stack[nesting].insertion || modifications_sequence_stack[nesting].deletion));
    if (modifications_sequence_stack[nesting].insertion || modifications_sequence_stack[nesting].deletion) {
        return;
    }

    if (nesting > 0 && (modifications_sequence_stack[nesting - 1].insertion || modifications_sequence_stack[nesting - 1].deletion)) {
        return;
    }

    --max_current_indels;

    bool was_insertion = false, was_deletion = false;
    for (int i = static_cast<int>(nesting) - 1; i >= std::max(static_cast<int>(nesting) - MIN_INS_DEL_DIST, 0); --i) {
        if (modifications_sequence_stack[i].insertion) {
            was_insertion = true;
        }
        else if (modifications_sequence_stack[i].deletion) {
            was_deletion = true;
        }
    }


    // try insertion

    if (!was_deletion) {
        // revert potential substitution correction
        modifications_sequence_stack[nesting].modification = sequence_modified[index_kmer];
        modifications_sequence_stack[nesting].insertion = true;

        if (index_kmer == 0) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                modifications_sequence_stack[nesting].insertion = false;
                return;
            }
            create_modification_path_with_indels_towards_5_prime(candidate_path_vector, nesting);
        }
        else {
            if (checked_changes == 0) {
                modifications_sequence_stack[nesting].insertion = false;
                return;
            }
            --checked_changes;
            std::string insertion_kmer = kmer.substr(1, kmer_length - 1);
            insertion_kmer += '0'; // dummy symbol, will be removed by extend_a_kmer function

            extend_a_kmer_5_prime_end(
                insertion_kmer,
                index_kmer - 1,
                candidate_path_vector,
                nesting + 1,
                checked_changes,
                max_remaining_changes,
                max_remaining_non_solid
                );
        }

        modifications_sequence_stack[nesting].insertion = false;
    }

    // try deletion

    if (!was_insertion) {
        std::string deletion_kmer(kmer);
        modifications_sequence_stack[nesting].deletion = true;
        // try different symbols for deletion
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                modifications_sequence_stack[nesting].deletion = false;
                return;
            }
            --checked_changes;

            deletion_kmer.back() = NEOCLEOTIDE[it_alter];

            float kmer_quality;
            if (query.query_text(deletion_kmer, kmer_quality)) {
                modifications_sequence_stack[nesting].quality = kmer_quality;
                modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                extend_a_kmer_5_prime_end(
                    deletion_kmer,
                    index_kmer, // next k-mer and get again the current one
                    candidate_path_vector,
                    nesting + 1,
                    checked_changes,
                    max_remaining_changes,
                    max_remaining_non_solid
                    );
            }
        }
        modifications_sequence_stack[nesting].deletion = false;
    }

    ++max_current_indels;
}



//----------------------------------------------------------------------
// Corrects indel situated on a 3' end of the read and traces recursively.
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::correct_indel_3_prime(const std::string & kmer, std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, std::size_t max_remaining_non_solid) {
    if (!modifications_sequence_stack[nesting].potential_indel) {
        return;
    }

    modifications_sequence_stack[nesting].potential_indel = false;
    --currently_proposed_indels;

    if (max_current_indels == 0) {
        return;
    }

    assert(!(modifications_sequence_stack[nesting].insertion || modifications_sequence_stack[nesting].deletion));
    if (modifications_sequence_stack[nesting].insertion || modifications_sequence_stack[nesting].deletion) {
        return;
    }

    if (nesting > 0 && (modifications_sequence_stack[nesting - 1].insertion || modifications_sequence_stack[nesting - 1].deletion)) {
        return;
    }

    --max_current_indels;

    bool was_insertion = false, was_deletion = false;
    for (int i = static_cast<int>(nesting) - 1; i >= std::max(static_cast<int>(nesting) - MIN_INS_DEL_DIST, 0); --i) {
        if (modifications_sequence_stack[i].insertion) {
            was_insertion = true;
        }
        else if (modifications_sequence_stack[i].deletion) {
            was_deletion = true;
        }
    }


    // try insertion

    if (!was_deletion) {
        // revert potential substitution correction
        modifications_sequence_stack[nesting].modification = sequence_modified[index_kmer + kmer_length - 1];
        modifications_sequence_stack[nesting].insertion = true;

        if (index_kmer == read_length - kmer_length) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                modifications_sequence_stack[nesting].insertion = false;
                return;
            }
            create_modification_path_with_indels_towards_3_prime(candidate_path_vector, nesting);
        }
        else {
            if (checked_changes == 0) {
                modifications_sequence_stack[nesting].insertion = false;
                return;
            }
            --checked_changes;
            std::string insertion_kmer("0"); // dummy symbol, will be removed by extend_a_kmer function
            insertion_kmer += kmer.substr(0, kmer_length - 1);

            extend_a_kmer_3_prime_end(
                insertion_kmer,
                index_kmer + 1,
                candidate_path_vector,
                nesting + 1,
                checked_changes,
                max_remaining_changes,
                max_remaining_non_solid
                );
        }

        modifications_sequence_stack[nesting].insertion = false;
    }

    // try deletion

    if (!was_insertion) {
        std::string deletion_kmer(kmer);
        modifications_sequence_stack[nesting].deletion = true;
        // try different symbols for deletion
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                modifications_sequence_stack[nesting].deletion = false;
                return;
            }
            --checked_changes;

            deletion_kmer.back() = NEOCLEOTIDE[it_alter];

            float kmer_quality;
            if (query.query_text(deletion_kmer, kmer_quality)) {
                modifications_sequence_stack[nesting].quality = kmer_quality;
                modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                extend_a_kmer_3_prime_end(
                    deletion_kmer,
                    index_kmer, // next k-mer and get again the current one
                    candidate_path_vector,
                    nesting + 1,
                    checked_changes,
                    max_remaining_changes,
                    max_remaining_non_solid
                    );
            }
        }
        modifications_sequence_stack[nesting].deletion = false;
    }

    ++max_current_indels;
}



//----------------------------------------------------------------------
// Proceeds correction towards 5' end.
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid) {
#ifdef LIMIT_MODIFICATIONS
    if (max_remaining_changes == -1) {
        return;
    }
#endif

    // generate a new k-mer
    std::string kmer_new(kmer.substr(0, kmer_length - 1));
    kmer_new = sequence_modified[index_kmer] + kmer_new;

    const bool is_low_quality_base = is_low_LUT[quality_score[index_kmer]];

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (!is_low_quality_base && query.query_text(kmer_new, kmer_quality) == true) {
        modifications_sequence_stack[nesting].quality = kmer_quality;
        modifications_sequence_stack[nesting].modification = kmer_new[0];

        // if this k-mer is the first k-mer in a read
        // running extend_a_kmer_5_prime_end is not needed any more
        if (index_kmer == 0) {
            if (candidate_path_vector.size() > MAX_EXTEND_CORRECTION_PATHS) {
                return;
            }
            create_modification_path_with_indels_towards_5_prime(candidate_path_vector, nesting);
        }
        else {
            extend_a_kmer_5_prime_end(
                    kmer_new,
                    index_kmer - 1,
                    candidate_path_vector,
                    nesting + 1,
                    checked_changes,
                    max_remaining_changes,
                    max_remaining_non_solid
                    );
        }
    }
    else {
        bool solid_found = false;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                return;
            }
            --checked_changes;
            // not equal to the original character
            if (sequence_modified[index_kmer] != NEOCLEOTIDE[it_alter] || is_low_quality_base) {
                // make a change
                kmer_new[0] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query.query_text(kmer_new, kmer_quality) == true) {
                    solid_found = true;

                    modifications_sequence_stack[nesting].quality = kmer_quality;
                    modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the first k-mer in a read
                    // running extend_a_kmer_5_prime_end is not needed any more
                    if (index_kmer == 0) {
                        if (candidate_path_vector.size() > MAX_EXTEND_CORRECTION_PATHS) {
                            return;
                        }
                        if (!create_modification_path_with_indels_towards_5_prime(candidate_path_vector, nesting)) {
                            check_is_new_potential_indel_5_prime(nesting, index_kmer, true);
                        }
                    }
                    else {

                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_5_prime_end(
                            kmer_new,
                            index_kmer - 1,
                            candidate_path_vector,
                            nesting + 1,
                            checked_changes,
                            (sequence_modified[index_kmer] != NEOCLEOTIDE[it_alter] ? max_remaining_changes - 1 : max_remaining_changes),
                            max_remaining_non_solid
                            );
                        check_is_new_potential_indel_5_prime(nesting, index_kmer);
                    }
                    correct_indel_5_prime(kmer_new, index_kmer, candidate_path_vector, nesting, checked_changes, max_remaining_changes, max_remaining_non_solid);
                }
            }
        }

        // try to accept non-solid k-mer
        if (!solid_found) {
            if (max_remaining_non_solid > 0) {
                kmer_new.front() = sequence_modified[index_kmer];

                modifications_sequence_stack[nesting].quality = 0;
                modifications_sequence_stack[nesting].modification = kmer_new[0];

                // if this k-mer is the first k-mer in a read
                // running extend_a_kmer_5_prime_end is not needed any more
                if (index_kmer == 0) {
                    if (candidate_path_vector.size() > MAX_EXTEND_CORRECTION_PATHS) {
                        return;
                    }
                    if (!create_modification_path_with_indels_towards_5_prime(candidate_path_vector, nesting)) {
                        check_is_new_potential_indel_5_prime(nesting, index_kmer);
                    }
                }
                else {
                    extend_a_kmer_5_prime_end(
                        kmer_new,
                        index_kmer - 1,
                        candidate_path_vector,
                        nesting + 1,
                        checked_changes,
                        max_remaining_changes,
                        max_remaining_non_solid - 1
                        );
                }
                correct_indel_5_prime(kmer_new, index_kmer, candidate_path_vector, nesting, checked_changes, max_remaining_changes, max_remaining_non_solid);
            }
        }
    }
}



template<>
inline void C_correct_read<false>::extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid) {
#ifdef LIMIT_MODIFICATIONS
    if (max_remaining_changes == -1) {
        return;
    }
#endif

    // generate a new k-mer
    std::string kmer_new(kmer.substr(0, kmer_length - 1));
    kmer_new = sequence_modified[index_kmer] + kmer_new;

    const bool is_low_quality_base = is_low_LUT[quality_score[index_kmer]];

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (!is_low_quality_base && query.query_text(kmer_new, kmer_quality) == true) {
        modifications_sequence_stack[nesting].quality = kmer_quality;
        modifications_sequence_stack[nesting].modification = kmer_new[0];

        // if this k-mer is the first k-mer in a read
        // running extend_a_kmer_5_prime_end is not needed any more
        if (index_kmer == 0) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                return;
            }
            create_modification_path_towards_5_prime(candidate_path_vector, nesting);
        }
        else {
            extend_a_kmer_5_prime_end(
                kmer_new,
                index_kmer - 1,
                candidate_path_vector,
                nesting + 1,
                checked_changes,
                max_remaining_changes,
                max_remaining_non_solid
                );
        }
    }
    else {
        bool solid_found = false;
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                return;
            }
            --checked_changes;

            // not equal to the original character
            if (sequence_modified[index_kmer] != NEOCLEOTIDE[it_alter] || is_low_quality_base) {
                // make a change
                kmer_new[0] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query.query_text(kmer_new, kmer_quality) == true) {
                    solid_found = true;
                    modifications_sequence_stack[nesting].quality = kmer_quality;
                    modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the first k-mer in a read
                    // running extend_a_kmer_5_prime_end is not needed any more
                    if (index_kmer == 0) {
                        if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                            return;
                        }
                        create_modification_path_towards_5_prime(candidate_path_vector, nesting);
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_5_prime_end(
                            kmer_new,
                            index_kmer - 1,
                            candidate_path_vector,
                            nesting + 1,
                            checked_changes,
                            (sequence_modified[index_kmer] != NEOCLEOTIDE[it_alter] ? max_remaining_changes - 1 : max_remaining_changes),
                            max_remaining_non_solid
                            );
                    }
                }
            }
        }

        // try to accept non-solid k-mer
        if (!solid_found) {
            if (max_remaining_non_solid > 0) {
                kmer_new.front() = sequence_modified[index_kmer];

                modifications_sequence_stack[nesting].quality = 0;
                modifications_sequence_stack[nesting].modification = kmer_new[0];

                // if this k-mer is the first k-mer in a read
                // running extend_a_kmer_5_prime_end is not needed any more
                if (index_kmer == 0) {
                    if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                        return;
                    }
                    create_modification_path_towards_5_prime(candidate_path_vector, nesting);
                }
                else {
                    extend_a_kmer_5_prime_end(
                        kmer_new,
                        index_kmer - 1,
                        candidate_path_vector,
                        nesting + 1,
                        checked_changes,
                        max_remaining_changes,
                        max_remaining_non_solid - 1
                        );
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Proceeds correction towards 3' end.
//----------------------------------------------------------------------

template<>
inline void C_correct_read<true>::extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid) {
#ifdef LIMIT_MODIFICATIONS
    if (max_remaining_changes == -1) {
        return;
    }
#endif

    // generate a new k-mer
    std::string kmer_new(kmer.substr(1, kmer_length - 1));
    kmer_new.push_back(sequence_modified[index_kmer + kmer_length - 1]);

    const bool is_low_quality_base = is_low_LUT[quality_score[index_kmer + kmer_length - 1]];

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (!is_low_quality_base && query.query_text(kmer_new, kmer_quality) == true) {
        modifications_sequence_stack[nesting].quality = kmer_quality;
        modifications_sequence_stack[nesting].modification = kmer_new[kmer_length - 1];

        // if this k-mer is the last k-mer in a read
        // running extend_a_kmer_3_prime_end is not needed any more
        if (index_kmer == (read_length - kmer_length)) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                return;
            }
            create_modification_path_with_indels_towards_3_prime(candidate_path_vector, nesting);
        }
        else {
            extend_a_kmer_3_prime_end(
                kmer_new,
                index_kmer + 1,
                candidate_path_vector,
                nesting + 1,
                checked_changes,
                max_remaining_changes,
                max_remaining_non_solid
                );
        }
    }
    else {
        bool solid_found = false;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                return;
            }
            --checked_changes;

            // not equal to the original character
            if (sequence_modified[index_kmer + kmer_length - 1] != NEOCLEOTIDE[it_alter] || is_low_quality_base) {
                // make a change
                kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query.query_text(kmer_new, kmer_quality) == true) {
                    solid_found = true;

                    modifications_sequence_stack[nesting].quality = kmer_quality;
                    modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the last k-mer in a read
                    // running extend_a_kmer_3_prime_end is not needed any more
                    if (index_kmer == (read_length - kmer_length)) {
                        if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                            return;
                        }
                        if (!create_modification_path_with_indels_towards_3_prime(candidate_path_vector, nesting)) {
                            check_is_new_potential_indel(nesting, index_kmer + kmer_length - 1, true);
                        }
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_3_prime_end(
                            kmer_new,
                            index_kmer + 1,
                            candidate_path_vector,
                            nesting + 1,
                            checked_changes,
                            sequence_modified[index_kmer + kmer_length - 1] != NEOCLEOTIDE[it_alter] ? max_remaining_changes - 1 : max_remaining_changes,
                            max_remaining_non_solid
                            );
                        check_is_new_potential_indel(nesting, index_kmer + kmer_length - 1);
                    }
                    correct_indel_3_prime(kmer_new, index_kmer, candidate_path_vector, nesting, checked_changes, max_remaining_changes, max_remaining_non_solid);
                }
            }
        }

        // try to accept non-solid k-mer
        if(!solid_found) {
            if (max_remaining_non_solid > 0) {
                kmer_new.back() = sequence_modified[index_kmer + kmer_length - 1];

                modifications_sequence_stack[nesting].quality = 0;
                modifications_sequence_stack[nesting].modification = kmer_new[kmer_length - 1];

                // if this k-mer is the last k-mer in a read
                // running extend_a_kmer_3_prime_end is not needed any more
                if (index_kmer == (read_length - kmer_length)) {
                    if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                        return;
                    }
                    if (!create_modification_path_with_indels_towards_3_prime(candidate_path_vector, nesting)) {
                        check_is_new_potential_indel(nesting, index_kmer + kmer_length - 1);
                    }
                }
                else {
                    extend_a_kmer_3_prime_end(
                        kmer_new,
                        index_kmer + 1,
                        candidate_path_vector,
                        nesting + 1,
                        checked_changes,
                        max_remaining_changes,
                        max_remaining_non_solid - 1
                        );
                }
                correct_indel_3_prime(kmer_new, index_kmer, candidate_path_vector, nesting, checked_changes, max_remaining_changes, max_remaining_non_solid);
            }
        }
    }
}



template<>
inline void C_correct_read<false>::extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid) {
#ifdef LIMIT_MODIFICATIONS
    if (max_remaining_changes == -1) {
        return;
    }
#endif

    // generate a new k-mer
    std::string kmer_new(kmer.substr(1, kmer_length - 1));
    kmer_new = kmer_new + sequence_modified[index_kmer + kmer_length - 1];

    const bool is_low_quality_base = is_low_LUT[quality_score[index_kmer + kmer_length - 1]];

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (!is_low_quality_base && query.query_text(kmer_new, kmer_quality) == true) {
        modifications_sequence_stack[nesting].quality = kmer_quality;
        modifications_sequence_stack[nesting].modification = kmer_new[kmer_length - 1];

        // if this k-mer is the last k-mer in a read
        // running extend_a_kmer_3_prime_end is not needed any more
        if (index_kmer == (read_length - kmer_length)) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                return;
            }
            create_modification_path_towards_3_prime(candidate_path_vector, index_kmer, nesting);
        }
        else {
            extend_a_kmer_3_prime_end(
                kmer_new,
                index_kmer + 1,
                candidate_path_vector,
                nesting + 1,
                checked_changes,
                max_remaining_changes,
                max_remaining_non_solid
                );
        }
    }
    else {
        bool solid_found = false;
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            if (checked_changes == 0) {
                return;
            }
            --checked_changes;

            // not equal to the original character
            if (sequence_modified[index_kmer + kmer_length - 1] != NEOCLEOTIDE[it_alter] || is_low_quality_base) {
                // make a change
                kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query.query_text(kmer_new, kmer_quality) == true) {
                    solid_found = true;
                    modifications_sequence_stack[nesting].quality = kmer_quality;
                    modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the last k-mer in a read
                    // running extend_a_kmer_3_prime_end is not needed any more
                    if (index_kmer == (read_length - kmer_length)) {
                        if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                            return;
                        }
                        create_modification_path_towards_3_prime(candidate_path_vector, index_kmer, nesting);
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_3_prime_end(
                            kmer_new,
                            index_kmer + 1,
                            candidate_path_vector,
                            nesting + 1,
                            checked_changes,
                            sequence_modified[index_kmer + kmer_length - 1] != NEOCLEOTIDE[it_alter] ? max_remaining_changes - 1 : max_remaining_changes,
                            max_remaining_non_solid
                            );
                    }
                }
            }
        }

        // try to accept non-solid k-mer
        if (!solid_found) {
            if (max_remaining_non_solid > 0) {
                kmer_new.back() = sequence_modified[index_kmer + kmer_length - 1];

                modifications_sequence_stack[nesting].quality = 0;
                modifications_sequence_stack[nesting].modification = kmer_new[kmer_length - 1];

                // if this k-mer is the last k-mer in a read
                // running extend_a_kmer_3_prime_end is not needed any more
                if (index_kmer == (read_length - kmer_length)) {
                    if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                        return;
                    }
                    create_modification_path_towards_3_prime(candidate_path_vector, index_kmer, nesting);
                }
                else {
                    extend_a_kmer_3_prime_end(
                        kmer_new,
                        index_kmer + 1,
                        candidate_path_vector,
                        nesting + 1,
                        checked_changes,
                        max_remaining_changes,
                        max_remaining_non_solid - 1
                        );
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Extends read towards 5'.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline void C_correct_read<CORRECT_INDEL>::extend_out_left(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, std::size_t max_remaining_non_solid, bool& extension_success) {
    // generate a new k-mer
    std::string kmer_new(kmer.substr(0, kmer_length - 1));
    kmer_new = '0' + kmer_new;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // generate kmer_new
        kmer_new[0] = NEOCLEOTIDE[it_alter];

        // kmer_new is solid
        float kmer_quality;
        bool solid_kmer = query.query_text(kmer_new, kmer_quality);
        if (solid_kmer || max_remaining_non_solid > 0) {
            // if current num_extend = extend_amount
            // running extend_out_left is not needed any more
            if ((num_extend + 1) == extend_amount) {
                extension_success = true;
                break;
            }
            else {
                // trace this kmer recursively
                extend_out_left(
                    kmer_new,
                    num_extend + 1,
                    extend_amount,
                    solid_kmer ? max_remaining_non_solid : max_remaining_non_solid - 1,
                    extension_success
                    );
                if (extension_success) {
                    break;
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Extends read towards 3'.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline void C_correct_read<CORRECT_INDEL>::extend_out_right(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, std::size_t max_remaining_non_solid, bool& extension_success) {
    // generate a new k-mer
    std::string kmer_new(kmer.substr(1, kmer_length - 1));
    kmer_new = kmer_new + '0';

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // generate kmer_new
        kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_new is solid
        float kmer_quality;
        bool solid_kmer = query.query_text(kmer_new, kmer_quality);
        if (solid_kmer || max_remaining_non_solid > 0) {
            // if current num_extend = extend_amount
            // running extend_out_right is not needed any more
            if ((num_extend + 1) == extend_amount) {
                extension_success = true;
                break;
            }
            else {
                // trace this kmer recursively
                extend_out_right(
                    kmer_new,
                    num_extend + 1,
                    extend_amount,
                    solid_kmer ? max_remaining_non_solid : max_remaining_non_solid - 1,
                    extension_success
                    );
                if (extension_success) {
                    break;
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Converts the quality indicator into an error probability.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
double C_correct_read<CORRECT_INDEL>::get_symbol_probability(std::size_t pos) {
    return quality_to_probability_LUT[quality_score[pos]];
}



//----------------------------------------------------------------------
// Creates vector of modifications basing on modifications buffer reading it towards 5'.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::create_modification_path_towards_5_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting) {
    // apply modifications
    C_candidate_path candidate_path;
    double kmers_quality = 0.0;
    for (std::size_t it = 0; it <= correction_nesting; ++it) {
        if (modifications_sequence_stack[it].modification != sequence_modified[correction_nesting - it]) {
            std::size_t modification_pos = correction_nesting - it;
            candidate_path.add_substitution(modification_pos, modifications_sequence_stack[it].modification, get_symbol_probability(modification_pos));
        }
        kmers_quality += modifications_sequence_stack[it].quality;
    }
    candidate_path.rate(kmers_quality, (correction_nesting + 1));

    return perform_extend_out_left(candidate_path, candidate_path_vector);
}



template<bool CORRECT_INDEL>
bool C_correct_read<CORRECT_INDEL>::create_modification_path_with_indels_towards_5_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting) {
    std::size_t current_modified_pos = current_region_index_start; // corrected symbol position in the input and output region; as std::string insert and erase method moves symbols on the right of the symbol, the only one index is required

    // apply modifications
    C_candidate_path candidate_path;
    double kmers_quality = 0.0;
    for (std::size_t it = 0; it <= correction_nesting; ++it) {
        assert(!(modifications_sequence_stack[it].insertion && modifications_sequence_stack[it].deletion));

        if (modifications_sequence_stack[it].insertion) {
            candidate_path.add_insertion(current_modified_pos);
        }

        else if (modifications_sequence_stack[it].deletion) {
            candidate_path.add_deletion(current_modified_pos + 1, modifications_sequence_stack[it].modification); // in std::string insert method ads symbol after n-th position, so add 1
        }
        else if (modifications_sequence_stack[it].modification != sequence_modified[current_modified_pos]) {
            candidate_path.add_substitution(current_modified_pos, modifications_sequence_stack[it].modification, get_symbol_probability(current_modified_pos));
        }

        if (!modifications_sequence_stack[it].insertion) {
            kmers_quality += modifications_sequence_stack[it].quality;
        }
        assert((it == correction_nesting && current_modified_pos == 0) || it < correction_nesting);
        if (!modifications_sequence_stack[it].deletion) {
            --current_modified_pos;
        }
    }

    candidate_path.rate(kmers_quality, (correction_nesting + 1 - candidate_path.get_read_length_change()));

    return perform_extend_out_left(candidate_path, candidate_path_vector);
}



//----------------------------------------------------------------------
// Creates vector of modifications basing on modifications buffer reading it towards 3'.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::create_modification_path_towards_3_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t current_kmer_pos, std::size_t correction_nesting) {
    const std::size_t first_symbol_pos = current_kmer_pos + kmer_length - 1 - correction_nesting;

    // apply modifications
    C_candidate_path candidate_path;
    double kmers_quality = 0.0;
    for (std::size_t it = 0; it <= correction_nesting; ++it) {
        std::size_t modification_pos = first_symbol_pos + it;
        if (modifications_sequence_stack[it].modification != sequence_modified[first_symbol_pos + it]) {
            candidate_path.add_substitution(modification_pos, modifications_sequence_stack[it].modification, get_symbol_probability(modification_pos));
        }
        kmers_quality += modifications_sequence_stack[it].quality;
    }
    candidate_path.rate(kmers_quality, (correction_nesting + 1));

    return perform_extend_out_right(candidate_path, candidate_path_vector);
}



template<bool CORRECT_INDEL>
bool C_correct_read<CORRECT_INDEL>::create_modification_path_with_indels_towards_3_prime(std::vector<C_candidate_path>& candidate_path_vector, std::size_t correction_nesting) {
    std::size_t current_modified_output_pos = current_region_index_start + kmer_length - 1; // corrected symbol position in the output region
    std::size_t current_modified_input_pos = current_region_index_start + kmer_length - 1; // corrected symbol position in the input region

    // apply modifications
    C_candidate_path candidate_path;
    double kmers_quality = 0.0;
    for (std::size_t it = 0; it <= correction_nesting; ++it) {
        assert(!(modifications_sequence_stack[it].insertion && modifications_sequence_stack[it].deletion));

        if (modifications_sequence_stack[it].insertion) {
            candidate_path.add_insertion(current_modified_output_pos);
        }
        else if (modifications_sequence_stack[it].deletion) {
            candidate_path.add_deletion(current_modified_output_pos, modifications_sequence_stack[it].modification);
        }
        else if (modifications_sequence_stack[it].modification != sequence_modified[current_modified_input_pos]) {
            candidate_path.add_substitution(current_modified_output_pos, modifications_sequence_stack[it].modification, get_symbol_probability(current_modified_input_pos));
        }

        if (!modifications_sequence_stack[it].insertion) {
            kmers_quality += modifications_sequence_stack[it].quality;

            ++current_modified_output_pos;
        }
        if (!modifications_sequence_stack[it].deletion) {
            ++current_modified_input_pos;
        }
    }

    candidate_path.rate(kmers_quality, (correction_nesting + 1 - candidate_path.get_read_length_change()));

    return perform_extend_out_right(candidate_path, candidate_path_vector);
}




//----------------------------------------------------------------------
// Creates vector of modifications basing on modifications buffer reading it towards 3'.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::create_modification_path_towards_3_prime_internal(std::vector<C_candidate_path>& candidate_path_vector, std::size_t current_kmer_pos, std::size_t correction_nesting) {
    const std::size_t first_symbol_pos = current_kmer_pos + kmer_length - 1 - correction_nesting;

    // apply modifications
    C_candidate_path candidate_path;
    double kmers_quality = 0.0;
    for (std::size_t it = 0; it <= correction_nesting; ++it) {
        std::size_t modification_pos = first_symbol_pos + it;
        if (modifications_sequence_stack[it].modification != sequence_modified[first_symbol_pos + it]) {
            candidate_path.add_substitution(modification_pos, modifications_sequence_stack[it].modification, get_symbol_probability(modification_pos));
        }
        kmers_quality += modifications_sequence_stack[it].quality;
    }
    candidate_path.rate(kmers_quality, (correction_nesting + 1));

    return perform_extend_right_inner_region(candidate_path, candidate_path_vector, current_kmer_pos);
}



//----------------------------------------------------------------------
// Creates vector of modifications basing on modifications buffer reading it towards 3'.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::create_modification_path_with_indels_towards_3_prime_internal(std::vector<C_candidate_path>& candidate_path_vector, std::size_t current_kmer_pos, std::size_t correction_nesting) {
    std::size_t current_modified_output_pos = current_region_index_start + kmer_length - 1; // corrected symbol position in the output region
    std::size_t current_modified_input_pos = current_region_index_start + kmer_length - 1; // corrected symbol position in the input region

    // apply modifications
    C_candidate_path candidate_path;
    double kmers_quality = 0.0;
    for (std::size_t it = 0; it <= correction_nesting; ++it) {
        assert(!(modifications_sequence_stack[it].insertion && modifications_sequence_stack[it].deletion));

        if (modifications_sequence_stack[it].insertion) {
            candidate_path.add_insertion(current_modified_output_pos);
        }
        else if (modifications_sequence_stack[it].deletion) {
            candidate_path.add_deletion(current_modified_output_pos, modifications_sequence_stack[it].modification);
        }
        else if (modifications_sequence_stack[it].modification != sequence_modified[current_modified_input_pos]) {
            candidate_path.add_substitution(current_modified_output_pos, modifications_sequence_stack[it].modification, get_symbol_probability(current_modified_input_pos));
        }

        if (!modifications_sequence_stack[it].insertion) {
            kmers_quality += modifications_sequence_stack[it].quality;

            ++current_modified_output_pos;
        }
        if (!modifications_sequence_stack[it].deletion) {
            ++current_modified_input_pos;
        }
    }

    candidate_path.rate(kmers_quality, (correction_nesting + 1 - candidate_path.get_read_length_change()));

    return perform_extend_right_inner_region(candidate_path, candidate_path_vector, current_kmer_pos);
}



//----------------------------------------------------------------------
// Creates vector of a single k-mer modifications basing on modifications buffer.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::create_modification_path_towards_3_prime_single_kmer(std::vector<C_candidate_path>& candidate_path_vector, const std::vector<std::size_t>& candidates_indexes, float kmer_quality) {
    C_candidate_path candidate_path;
    for (std::size_t it = 0; it < candidates_indexes.size(); ++it) {
        const std::size_t current_pos = candidates_indexes[it];
        if (sequence_modified[current_pos] != modifications_sequence_stack[it].modification) {
            candidate_path.add_substitution(current_pos, modifications_sequence_stack[it].modification, get_symbol_probability(current_pos));
        }
        candidate_path.add_substitution(current_pos, modifications_sequence_stack[it].modification, get_symbol_probability(current_pos));
    }
    candidate_path.rate(kmer_quality, 1U);

    return perform_extend_out_first_kmer(candidate_path, candidate_path_vector);
}



//----------------------------------------------------------------------
// Creates modification of a single k-mer vector.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline bool C_correct_read<CORRECT_INDEL>::create_modification_path_with_single_substitution(std::vector<C_candidate_path>& candidate_path_vector, float kmer_quality, std::size_t pos, char modification) {
    C_candidate_path candidate_path;

    candidate_path.add_substitution(pos, modification, get_symbol_probability(pos));
    candidate_path.rate(kmer_quality, 1U);

    return perform_extend_out_first_kmer(candidate_path, candidate_path_vector);
}



//----------------------------------------------------------------------
// Creates modification path with no modifications and does not extends such a path.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
inline void C_correct_read<CORRECT_INDEL>::create_modification_path_empty_no_extend(std::vector<C_candidate_path>& candidate_path_vector, float kmer_quality) {
    C_candidate_path candidate_path;
    candidate_path.rate(kmer_quality, 1ULL);
    candidate_path_vector.push_back(std::move(candidate_path));
}



//----------------------------------------------------------------------
// Applies changes into the read.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::modify_errors(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors) {
    std::vector<C_candidate_path>::iterator best_it_path = choose_best_correction(candidate_path_vector);
    if (best_it_path != candidate_path_vector.end()) {
        apply_path_to_read(*best_it_path, sequence_modified, quality_score, num_corrected_errors);
    }
}



//----------------------------------------------------------------------
// Applies changes into the read, but does not select it by rate.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::modify_errors_no_rate(C_candidate_path& candidate_path, std::size_t& num_corrected_errors) {
    apply_path_to_read(candidate_path, sequence_modified, quality_score, num_corrected_errors);
}



//----------------------------------------------------------------------
// Applies changes into the read with respect of modification position.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::modify_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2) {
    std::vector<C_candidate_path>::iterator best_it_path = choose_best_correction(candidate_path_vector);
    if (best_it_path != candidate_path_vector.end()) {
        apply_path_to_read(*best_it_path, sequence_modified, quality_score, num_corrected_errors1, num_corrected_errors2);
    }
}



//----------------------------------------------------------------------
// Applies changes into the read with corrected first k-mer, but does not select it by rate.
//----------------------------------------------------------------------

template<bool CORRECT_INDEL>
void C_correct_read<CORRECT_INDEL>::modify_errors_first_kmer_no_rate(C_candidate_path& candidate_path, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2) {
    apply_path_to_read(candidate_path, sequence_modified, quality_score, num_corrected_errors1, num_corrected_errors2);
}

#endif /* CORRECT_READ_HPP */
