/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#ifndef CORRECT_READ_HPP
#define	CORRECT_READ_HPP



#include "parse_args.hpp"
#include <kmc_api/kmc_file.h>
#include <kmc_api/kmer_api.h>
#include <vector>
#include <fstream>



//#define USE_KMER_MEDIAN



//----------------------------------------------------------------------
// C_candidate_path
//----------------------------------------------------------------------

class C_candidate_path {
public:
    typedef std::pair<std::size_t, char> Single_mod;
    // std::string last_kmer;
    // modified positions and modification results
    std::vector< Single_mod > modified_bases;

    double kmers_quality;

#ifdef USE_KMER_MEDIAN
    std::vector<float> covering_kmers_weight_vector;
#else
    double covering_kmers_weight;
#endif

    // constructors

    C_candidate_path() : kmers_quality(0.0)
#ifndef USE_KMER_MEDIAN
    , covering_kmers_weight(0.0)
#endif
    {
    };

    C_candidate_path(const float _covering_kmer_weight, const float _kmer_quality,
            const std::size_t n_modifications, const char modifications[], const std::vector<std::size_t> mod_pos) :
    kmers_quality(_kmer_quality),
    covering_kmers_weight(_covering_kmer_weight) {
        modified_bases.reserve(n_modifications);
        for (std::size_t it = 0; it < n_modifications; ++it) {
            Single_mod mod(mod_pos[it], modifications[it]);
            modified_bases.push_back(mod);
        }
    };

    void clear_path();
};

template <std::size_t MAX_MODIFICATIONS>
class C_fast_candidate_path {
public:
    char modifications[MAX_MODIFICATIONS];
};

class C_modification_with_quality {
public:
    float quality;
    char modification;
};

class C_correct_read {
public:
    // variables
    std::size_t read_length;
    std::size_t quality_score_offset;
    std::size_t kmer_length;
    std::size_t max_extension;
    const std::size_t max_read_length;

    // constructors

    C_correct_read(std::size_t _quality_score_offset, std::size_t _kmer_length, std::size_t _max_extension, const std::size_t _max_read_length,
            const std::string& _read_file_name,
            std::size_t& _num_corrected_errors_step1_1, std::size_t& _num_corrected_errors_step1_2,
            std::size_t& _num_corrected_errors_step1_3, std::size_t& _num_corrected_errors_step2_1,
            std::size_t& _num_corrected_errors_step2_2, CKMCFile& _kmc_file) :
    read_length(0),
    quality_score_offset(_quality_score_offset),
    kmer_length(_kmer_length),
    max_extension(_max_extension),
    max_read_length(_max_read_length),
    sequence_modification(_max_read_length, '0'),
    kmc_file(_kmc_file),
    kmer_api(_kmer_length),
    num_corrected_errors_step1_1(_num_corrected_errors_step1_1),
    num_corrected_errors_step1_2(_num_corrected_errors_step1_2),
    num_corrected_errors_step1_3(_num_corrected_errors_step1_3),
    num_corrected_errors_step2_1(_num_corrected_errors_step2_1),
    num_corrected_errors_step2_2(_num_corrected_errors_step2_2) {
    }

    void correct_errors_in_a_read_fastq();

    std::string header;
    std::string sequence;
    std::string connector;
    std::string quality_score;

    std::string sequence_modification;

private:
    typedef C_candidate_path::Single_mod Single_mod;

    CKMCFile& kmc_file;
    CKmerAPI kmer_api;

    std::size_t& num_corrected_errors_step1_1;
    std::size_t& num_corrected_errors_step1_2;
    std::size_t& num_corrected_errors_step1_3;
    std::size_t& num_corrected_errors_step2_1;
    std::size_t& num_corrected_errors_step2_2;

    std::string sequence_modified;

    // functions
    bool query_text(const std::string& kmer, float& kmer_quality);

    void correct_errors_between_solid_regions(const std::size_t index_start, const std::size_t index_end);
    void correct_errors_5_prime_end(const std::size_t index_start);
    void correct_errors_3_prime_end(const std::size_t index_start);
    void correct_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector);
    void extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector);
    void extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t max_remaining_changes, C_modification_with_quality modifications[], std::size_t nesting, std::size_t& checked_changes);
    void extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t max_remaining_changes, C_modification_with_quality modifications[], std::size_t nesting, std::size_t& checked_changes);
    void perform_extend_out_left(std::string& sequence_tmp, C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector_tmp_tmp);
    void perform_extend_out_right(std::string& sequence_tmp, C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector_tmp_tmp);
    void extend_out_left(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, bool& extension_success);
    void extend_out_right(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, bool& extension_success);
    void check_first_kmer(std::string& kmer, C_fast_candidate_path<MAX_CHECK_FIRST_KMER_NESTING>& candidate_path_in, const std::vector<std::size_t>& candidates_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t index);
    void solid_first_kmer(const C_candidate_path& candidate_path, bool& extension_success);
    void extend_first_kmer_to_right(C_candidate_path& candidate_path_in, bool& correction_success);
    double convert_quality_to_probability(char c);

    std::vector<C_candidate_path>::iterator choose_best_correction(std::vector<C_candidate_path>& candidate_path_vector);
    void modify_errors(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors);
    void modify_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2);

#ifdef USE_KMER_MEDIAN
    float median(std::vector<float>& kmer_qualities);
#endif
};
#endif	/* CORRECT_READ_HPP */

