/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#include "correct_read.hpp"
#include "Log.h"
#include <algorithm>



//----------------------------------------------------------------------
// Determines erroneous regions and calls their correction.
//----------------------------------------------------------------------

void C_correct_read::correct_errors_in_a_read_fastq() {
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
        std::string current_kmer(sequence.substr(it_kmer, kmer_length));

        // k-mer is solid
        float kmer_quality;
        if (query_text(current_kmer, kmer_quality) == true) {
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
                    std::cout << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
                    Log::get_stream() << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
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
    // STEP 0-1: remove short solid regions beside short non-solid regions
    //--------------------------------------------------
    /*
    if (solid_regions.size() > 0) {
       std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
       solid_regions_tmp.push_back(solid_regions[0]);

       if (solid_regions.size() > 1) {
          for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
             // short non-solid region to the left
             if ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < kmer_length) {
                if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) < MIN_SOLID_LENGTH) {
                   // (non-solid region length < kmer_length) &&
                   // (solid_regions[it_region] is too short) &&
                   // (solid_regions[it_region].second is not the last base of the read)
                   // -> remove solid_regions[it_region + 1]
                   // k = 10
                   // SSSSSNNNNNNNSNNSSSSS
                   //             ^ Remove it
                   // do nothing
                   if (solid_regions[it_region].second != (sequence.length() - 1)) {
                   }
                   else {
                      solid_regions_tmp.push_back(solid_regions[it_region]);
                   }
                }
                else {
                   solid_regions_tmp.push_back(solid_regions[it_region]);
                }
             }
             else {
                solid_regions_tmp.push_back(solid_regions[it_region]);
             }
          }
       }
       solid_regions = solid_regions_tmp;
    }
     */

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

        solid_regions = solid_regions_tmp;
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
        solid_regions = solid_regions_tmp;
    }

    //--------------------------------------------------
    // STEP 0-4: reduce the size of solid regions
    //--------------------------------------------------
    if (solid_regions.size() > 1) {
        for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
            // (length of a non-solid region < kmer_length) && (length of a non-solid region >= kmer_length - FP_SUSPECT_LENGTH(default: 1))
            if (((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < kmer_length) &&
                    ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) >= kmer_length - FP_SUSPECT_LENGTH)) {
                // length of the right solid region > FP_SUSPECT_LENGTH(default: 1)
                if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) > FP_SUSPECT_LENGTH) {
                    solid_regions[it_region].first += FP_SUSPECT_LENGTH;
                }

                // length of the left solid region > FP_SUSPECT_LENGTH(default: 1)
                if ((solid_regions[it_region - 1].second - solid_regions[it_region - 1].first + 1) > FP_SUSPECT_LENGTH) {
                    solid_regions[it_region - 1].second -= FP_SUSPECT_LENGTH;
                }
            }
        }
    }

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
                    if ((((std::size_t)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < QS_CUTOFF) ||
                            (((std::size_t)quality_score[it_adjust] - quality_score_offset) < QS_CUTOFF)
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
                    if ((((std::size_t)quality_score[it_adjust] - quality_score_offset) < QS_CUTOFF) ||
                            (((std::size_t)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < QS_CUTOFF)
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
                if (((std::size_t)quality_score[it_adjust] - quality_score_offset) < QS_CUTOFF) {
                    solid_regions[0].first = it_adjust + 1;
                }
            }
        }
    }

    //--------------------------------------------------
    // STEP 0-7: check whether a non-solid region < k still exists
    //--------------------------------------------------
    bool short_non_solid_region(false);
    if (solid_regions.size() > 1) {
        for (std::size_t it_sr = 1; it_sr < (solid_regions.size() - 1); it_sr++) {
            if ((solid_regions[it_sr].first - solid_regions[it_sr - 1].second) <= kmer_length) {
                short_non_solid_region = true;
                break;
            }
        }
    }

    //--------------------------------------------------
    // correct errors
    //--------------------------------------------------
    sequence_modified = sequence;

    if ((solid_regions.size() > 0) && (short_non_solid_region == false)) {
        //--------------------------------------------------
        // STEP 1-1: Correct errors between solid regions
        //--------------------------------------------------
        if (solid_regions.size() > 1) {
            // for each solid region
            for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
                if ((((solid_regions[it_region].first - 1) - (solid_regions[it_region - 1].second + 1)) + 1) >= kmer_length) {
                    correct_errors_between_solid_regions(
                            (solid_regions[it_region - 1].second + 1),
                            (solid_regions[it_region].first - 1)
                            );
                }
                else {
                }
            }
        }

        //--------------------------------------------------
        // STEP 1-2: Correct errors in the 5' end
        //--------------------------------------------------
        // number of solid regions is >= 1
        if (solid_regions.size() >= 1) {
            // the first solid region does not start from the 0-th k-mer in a read
            if (solid_regions[0].first > 0) {

                correct_errors_5_prime_end(solid_regions[0].first - 1);

            }
        }

        //--------------------------------------------------
        // STEP 1-3: Correct errors in the 3' end
        //--------------------------------------------------
        // number of solid regions is >= 1
        if (solid_regions.size() >= 1) {
            // the last solid region does not end in the last k-mer in a read
            if (solid_regions[solid_regions.size() - 1].second < (read_length - kmer_length)) {

                correct_errors_3_prime_end(solid_regions[solid_regions.size() - 1].second + 1);

            }
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
        std::vector<C_candidate_path> candidate_path_vector_tmp;

        correct_errors_first_kmer(candidate_path_vector_tmp);

        // candidiate_path_vector_tmp: differently modified versions of the first k-mer

        // filter some candidates by extending the first k-mer to the left
        std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

        if (candidate_path_vector_tmp.size() > 0) {
            // each path
            for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size(); it_candidates++) {
                // no modified path
                if (candidate_path_vector_tmp[it_candidates].modified_bases.size() == 0) {
                    candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
                }
                    // check the index of the first modified base
                    // extension is needed
                    //else if (candidate_path_vector_tmp[it_candidates].modified_bases[0].first < (MAX_EXTENSION - 1)) {
                else if (candidate_path_vector_tmp[it_candidates].modified_bases[0].first < (kmer_length - 1)) {
                    bool extension_success(false);
                    solid_first_kmer(
                            candidate_path_vector_tmp[it_candidates],
                            extension_success
                            );

                    if (extension_success == true) {
                        candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
                    }
                }
                    // extension is not needed
                else {
                    candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
                }
            }
        }

        // candidiate_path_vector_tmp_tmp: solid k-mers in candidate_path_vector_tmp

        // candidates in candidiate_path_vector_tmp_tmp are moved to candidate_path_vector_tmp
        candidate_path_vector_tmp = candidate_path_vector_tmp_tmp;
        candidate_path_vector_tmp_tmp.clear();

        //--------------------------------------------------
        // STEP 2-2: extend candidate paths to the right
        //--------------------------------------------------
        if (candidate_path_vector_tmp.size() > 0) {
            // each path
            for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size(); it_candidates++) {
                bool correction_success(false);

                extend_first_kmer_to_right(
                        candidate_path_vector_tmp[it_candidates],
                        correction_success
                        );

                // add this path to candidate_path_vector_tmp_tmp if its correction succeeds
                if (correction_success == true) {
                    candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
                }
            }
        }

        // candidiate_path_vector_tmp_tmp: successfully right extended candidates

        //--------------------------------------------------
        // STEP 2-3: choose a final one in candidate_path_vector_tmp_tmp if possible
        //--------------------------------------------------
        // compare quality scores of candidate paths
        // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1

        modify_errors_first_kmer(candidate_path_vector_tmp_tmp, num_corrected_errors_step2_1, num_corrected_errors_step2_2);
    }
}



//----------------------------------------------------------------------
// Corrects errors in a region situated between correct regions.
//----------------------------------------------------------------------

inline void C_correct_read::correct_errors_between_solid_regions(const std::size_t index_start, const std::size_t index_end) {
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
    // list of candidate paths
    std::vector<C_candidate_path> candidate_path_vector_tmp;

    // index of the k-mer that can be modified
    // k-mers that are overlapped with a solid regioin cannot be modified
    std::size_t index_last_mod(index_end - kmer_length + 1);

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query_text(kmer_initial, kmer_quality) == true) {
            // generate a new path
            C_candidate_path candidate_path;

            if (sequence_modified[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
                Single_mod pair_tmp(index_start + kmer_length - 1, NEOCLEOTIDE[it_alter]);

                candidate_path.modified_bases.push_back(pair_tmp);
            }
            candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
            candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
            candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

            // if this k-mer is the last k-mer that can be modified
            // running extend_a_kmer_right is not needed any more
            if (index_start == index_last_mod) {
                candidate_path_vector_tmp.push_back(candidate_path);
            }
            else {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_a_kmer(
                        kmer_initial,
                        index_start,
                        index_last_mod,
                        candidate_path,
                        candidate_path_vector_tmp
                        );
            }
        }
    }

    std::vector<C_candidate_path> candidate_path_vector;

    // check the solidness of k-mers between index_last_mod and index_end
    bool all_solid_wo_modification(false);

    // each candidate path
    for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); ++it_path) {
        if ((*it_path).modified_bases.size() == 0) {
            all_solid_wo_modification = true;
            break;
        }
        else {
            // checking is needed
            std::size_t index_last_modified_base((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

            if (index_last_modified_base > index_last_mod) {
                // generate a temporary sequence
                std::string sequence_tmp(sequence_modified);
                for (std::size_t it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
                    sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
                }

                // check k-mers
                std::size_t num_success(0);
                for (std::size_t it_check = index_last_mod; it_check <= index_last_modified_base; it_check++) {
                    std::string kmer_current(sequence_tmp.substr(it_check, kmer_length));

                    float kmer_quality;
                    if (query_text(kmer_current, kmer_quality) == true) {
                        num_success++;
                    }
                    else {
                        break;
                    }
                }

                if (num_success == (index_last_modified_base - index_last_mod + 1)) {
                    candidate_path_vector.push_back(*it_path);
                }
            }
                // checking is not needed
            else {
                candidate_path_vector.push_back(*it_path);
            }
        }
    }

    // all k-mers are solid without any modification
    // do nothing
    if (all_solid_wo_modification == true) {
    }
        // compare quality scores of candidate paths
        // if the number of paths in candidate_path_vector is larger than 1
    else {
        modify_errors(candidate_path_vector, num_corrected_errors_step1_1);
    }
}



//----------------------------------------------------------------------
// Corrects errors situated on a 5' end of the read
//----------------------------------------------------------------------

inline void C_correct_read::correct_errors_5_prime_end(const std::size_t index_start) {
    // |  non-solid  | 1st solid region
    // |--------------------------------------| read
    //         |-----|                          (index_start)-th k-mer
    //--------------------------------------------------
    // list of candidate paths
    std::vector<C_candidate_path> candidate_path_vector_tmp;

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[0] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query_text(kmer_initial, kmer_quality) == true) {
            // if this k-mer is the first k-mer in a read
            // running extend_a_kmer_5_prime_end is not needed any more
            if (index_start == 0) {
                // generate a new path
                C_candidate_path candidate_path;

                Single_mod pair_tmp(index_start, NEOCLEOTIDE[it_alter]);

                candidate_path.modified_bases.push_back(pair_tmp);

                candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
                candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
                candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif
                candidate_path_vector_tmp.push_back(candidate_path);
            }
            else if (index_start > 0) {
                const std::size_t max_remaining_changes = (NEOCLEOTIDE[it_alter] == sequence_modified[index_start]
                        ? (MAX_CHANGES_IN_REGION_RATIO * index_start) : (MAX_CHANGES_IN_REGION_RATIO * index_start) - 1);
                C_modification_with_quality modifications[max_read_length];
                modifications[0].modification = NEOCLEOTIDE[it_alter];
                modifications[0].quality = kmer_quality;
                // trace  this kmer recursively and update candidate_path_vector_tmp

                std::size_t checked_changes = 0;
                extend_a_kmer_5_prime_end(
                        kmer_initial,
                        index_start,
                        candidate_path_vector_tmp,
                        max_remaining_changes,
                        modifications,
                        1,
                        checked_changes
                        );
            }
        }
        else {
        }
    }

    std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

    // each candidate path
    for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); ++it_path) {
        std::string sequence_tmp(sequence_modified);
        perform_extend_out_left(sequence_tmp, *it_path, candidate_path_vector_tmp_tmp);
    }

    modify_errors(candidate_path_vector_tmp_tmp, num_corrected_errors_step1_2);
}



//----------------------------------------------------------------------
// Corrects errors situated on a 3' end of the read
//----------------------------------------------------------------------

inline void C_correct_read::correct_errors_3_prime_end(const std::size_t index_start) {
    //  last solid region | non-solid region |
    // --------------------------------------| read
    //               |-----|                   (index_start)-th k-mer
    //--------------------------------------------------
    // list of candidate paths
    std::vector<C_candidate_path> candidate_path_vector_tmp;

    // make an initial k-mer
    std::string kmer_initial(sequence_modified.substr(index_start, kmer_length));

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query_text(kmer_initial, kmer_quality) == true) {
            // if this k-mer is the last k-mer in a read
            // running extend_a_kmer_3_prime_end is not needed any more
            if (index_start == (read_length - kmer_length)) {
                // generate a new path
                C_candidate_path candidate_path;

                Single_mod pair_tmp(index_start + kmer_length - 1, NEOCLEOTIDE[it_alter]);

                candidate_path.modified_bases.push_back(pair_tmp);

                candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
                candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
                candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif
                candidate_path_vector_tmp.push_back(candidate_path);
            }
            else if (index_start < (read_length - kmer_length)) {
                const std::size_t max_remaining_changes = (NEOCLEOTIDE[it_alter] == sequence_modified[index_start]
                        ? (MAX_CHANGES_IN_REGION_RATIO * (read_length - index_start)) : (MAX_CHANGES_IN_REGION_RATIO * (read_length - index_start)) - 1);
                C_modification_with_quality modifications[max_read_length];
                modifications[0].modification = NEOCLEOTIDE[it_alter];
                modifications[0].quality = kmer_quality;
                // trace  this kmer recursively and update candidate_path_vector_tmp
                C_candidate_path temp_path;

                std::size_t checked_changes = 0;
                extend_a_kmer_3_prime_end(
                        kmer_initial,
                        index_start,
                        temp_path,
                        candidate_path_vector_tmp,
                        max_remaining_changes,
                        modifications,
                        1,
                        checked_changes
                        );
            }
        }
    }

    std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

    // each candidate path
    for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); ++it_path) {
        std::string sequence_tmp(sequence_modified);
        perform_extend_out_right(sequence_tmp, *it_path, candidate_path_vector_tmp_tmp);
    }

    modify_errors(candidate_path_vector_tmp_tmp, num_corrected_errors_step1_3);
}



//----------------------------------------------------------------------
// Corrects errors in the first k-mer of the read.
//----------------------------------------------------------------------

inline void C_correct_read::correct_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector) {
    std::string first_kmer(sequence_modified.substr(0, kmer_length));

    std::vector<std::size_t> low_qs_indexes;

    for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
        if (((std::size_t)quality_score[it_bases] - quality_score_offset) < QS_CUTOFF) {
            low_qs_indexes.push_back(it_bases);
        }
    }

    // correct errors if the number of low-quality bases is smaller than the threshold
    if ((low_qs_indexes.size() <= MAX_LOW_QS_BASES) && (low_qs_indexes.size() > 0)) {
        C_fast_candidate_path<MAX_CHECK_FIRST_KMER_NESTING> candidate_fast_path;

        std::string kmer(first_kmer);
        check_first_kmer(
                kmer,
                candidate_fast_path,
                low_qs_indexes,
                candidate_path_vector,
                0
                );

        // no candidate path is found
        if (candidate_path_vector.size() == 0) {
            for (std::size_t it_bases = 0; it_bases < kmer_length; it_bases++) {
                std::string kmer_tmp(first_kmer);

                // each alternative neocletide
                for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
                    // not equal to the original character
                    if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                        // generate a new k-mer
                        kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                        // add kmer_tmp to candidate_path_tmp if it is solid
                        float kmer_quality;
                        if (query_text(kmer_tmp, kmer_quality) == true) {
                            // generate a new candidate path
                            C_candidate_path candidate_path;

                            Single_mod pair_tmp(it_bases, NEOCLEOTIDE[it_alter]);
                            candidate_path.modified_bases.push_back(pair_tmp);

                            candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
                            candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
                            candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

                            candidate_path_vector.push_back(candidate_path);
                        }
                    }
                }
            }
        }
    }
        // no low-quality base or too many low-quality bases
    else {
        float kmer_quality;
        if (query_text(first_kmer, kmer_quality) == true) {
            C_candidate_path candidate_path;
            candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
            candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
            candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

            candidate_path_vector.push_back(candidate_path);
        }
        else {
            // (quality, position)
            std::vector<std::pair<char, std::size_t> > qualities_vector;
            qualities_vector.reserve(kmer_length);

            for (std::size_t it_qualities = 0; it_qualities < kmer_length; it_qualities++) {
                qualities_vector.push_back(std::pair<char, std::size_t>(quality_score[it_qualities], it_qualities));
            }

            std::sort(qualities_vector.begin(), qualities_vector.end());

            if (low_qs_indexes.size() == 0) {
                for (std::size_t i = 0; i < kmer_length; i++) {
                    std::size_t it_bases = qualities_vector[i].second;
                    std::string kmer_tmp(first_kmer);

                    // each alternative neocletide
                    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
                        // not equal to the original character
                        if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                            // generate a new k-mer
                            kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                            // add kmer_tmp to candidate_path_tmp if it is solid
                            float kmer_quality;
                            if (query_text(kmer_tmp, kmer_quality) == true) {
                                // generate a new candidate path
                                C_candidate_path candidate_path;

                                Single_mod pair_tmp(it_bases, NEOCLEOTIDE[it_alter]);
                                candidate_path.modified_bases.push_back(pair_tmp);

                                candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
                                candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
                                candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

                                candidate_path_vector.push_back(candidate_path);
                            }
                        }
                    }
                }
            }
            else {
                std::vector<std::size_t> candidates;

                for (std::size_t it = 0; it < std::min((std::size_t)low_qs_indexes.size(), (std::size_t)MAX_LOW_QS_INDEXES_COMB); ++it) {
                    candidates.push_back(qualities_vector[it].second);
                }

                if (candidates.size() > 0) {
                    C_fast_candidate_path<MAX_CHECK_FIRST_KMER_NESTING> candidate_path;

                    std::string kmer(first_kmer);
                    check_first_kmer(
                            kmer,
                            candidate_path,
                            candidates,
                            candidate_path_vector,
                            0
                            );
                }
            }
        }
    }

    if (candidate_path_vector.size() > MAX_FIRST_KMER_POSIBILITIES) {
        std::vector<C_candidate_path> candidate_path_vector_tmp_tmp(candidate_path_vector);

        candidate_path_vector.clear();

        // (probability, index)
        std::vector<std::pair<double, std::vector<C_candidate_path>::iterator> > probabilities;
        // each candidate path

        for(std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); ++it_path) {
            // each modification
#ifdef USE_KMER_MEDIAN
            double prob = it_path->covering_kmers_weight_vector[0];
#else
            double prob = it_path->covering_kmers_weight;
#endif

            for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
                // add quality scores of modified bases
                if (sequence_modified[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
                    prob *= convert_quality_to_probability(quality_score[(*it_path).modified_bases[it_mod].first]);
                }
            }

            probabilities.push_back(make_pair(prob, it_path));
        }

        std::sort(probabilities.begin(), probabilities.end());

        for (std::size_t it = 0; it < (std::size_t)MAX_FIRST_KMER_POSIBILITIES; ++it) {
            candidate_path_vector.push_back(*probabilities[probabilities.size() - it - 1].second);
        }
    }
}



//----------------------------------------------------------------------
// Tries to correct k-mer by changing one symbol.
//----------------------------------------------------------------------

inline void C_correct_read::check_first_kmer(std::string& kmer, C_fast_candidate_path<MAX_LOW_QS_INDEXES_COMB>& candidate_path_in, const std::vector<std::size_t>& candidates_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t index) {
    if (candidate_path_vector.size() >= MAX_FIRST_KMER_CORRECTION_PATHS) {
        return;
    }
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // make a new k-mer
        kmer[candidates_indexes[index]] = NEOCLEOTIDE[it_alter];

        candidate_path_in.modifications[index] = NEOCLEOTIDE[it_alter];

        if (index == candidates_indexes.size() - 1) {
            float kmer_quality;
            if (query_text(kmer, kmer_quality) == true) {
                C_candidate_path path_temp(COVERING_KMERS_WEIGHT, kmer_quality, candidates_indexes.size(), candidate_path_in.modifications, candidates_indexes);
                candidate_path_vector.push_back(path_temp);
            }
        }
        else {
            check_first_kmer(
                    kmer,
                    candidate_path_in,
                    candidates_indexes,
                    candidate_path_vector,
                    index + 1
                    );
        }
    }
}



//----------------------------------------------------------------------
// Checks if modified first k-mer can be extended.
//----------------------------------------------------------------------

inline void C_correct_read::solid_first_kmer(const C_candidate_path& candidate_path, bool& extension_success) {
    // index_smallest_modified
    std::size_t index_smallest_modified(candidate_path.modified_bases[0].first);

    // number of bases that should be extended
    std::size_t extend_amount;

    // applied the modified bases to first_kmer
    std::string first_kmer(sequence_modified.substr(0, kmer_length));
    for (std::size_t it_base = 0; it_base < candidate_path.modified_bases.size(); it_base++) {
        first_kmer[candidate_path.modified_bases[it_base].first] = candidate_path.modified_bases[it_base].second;
    }

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

    // generate an initial k-mer
    std::string kmer_initial(first_kmer.substr(0, kmer_length - 1));
    kmer_initial = '0' + kmer_initial;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // make a change
        kmer_initial[0] = NEOCLEOTIDE[it_alter];

        // kmer_initial is solid
        float kmer_quality;
        if (query_text(kmer_initial, kmer_quality) == true) {
            // if extend_amount == 1
            // running extend_out_left is not needed any more
            if (extend_amount == 1) {
                extension_success = true;
                break;
            }
            else if (!extension_success) {
                // trace  this kmer recursively and update candidate_path_vector_tmp
                extend_out_left(
                        kmer_initial,
                        1,
                        extend_amount,
                        extension_success
                        );
            }
        }
    }
}



//----------------------------------------------------------------------
// Proceeds correction from modified first k-mer towards 3' end.
//----------------------------------------------------------------------

inline void C_correct_read::extend_first_kmer_to_right(C_candidate_path& candidate_path_in, bool& correction_success) {
    // generate the first k-mer
    std::string first_kmer(sequence_modified.substr(0, kmer_length + 1));
    for (std::size_t it_base = 0; it_base < candidate_path_in.modified_bases.size(); it_base++) {
        first_kmer[candidate_path_in.modified_bases[it_base].first] = candidate_path_in.modified_bases[it_base].second;
    }

    // generate the second k-mer
    std::string second_kmer(first_kmer.substr(1, kmer_length));

    // list of candidate paths
    std::vector<C_candidate_path> candidate_path_vector_tmp;


    // second_kmer is solid
    float kmer_quality;
    if (query_text(second_kmer, kmer_quality) == true) {
        if ((read_length - kmer_length) == 1) {
            // if this k-mer is the last k-mer in a read
            // running extend_a_kmer_3_prime_end is not needed any more

            candidate_path_in.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
            candidate_path_in.covering_kmers_weight_vector.push_back(kmer_quality);
#else
            candidate_path_in.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

            candidate_path_vector_tmp.push_back(candidate_path_in);
        }
        else if ((read_length - kmer_length) > 1) {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            C_modification_with_quality modifications[max_read_length];
            modifications[0].modification = second_kmer[kmer_length - 1];
            modifications[0].quality = kmer_quality;

            std::size_t checked_changes = 0;
            extend_a_kmer_3_prime_end(
                    second_kmer,
                    1,
                    candidate_path_in,
                    candidate_path_vector_tmp,
                    (MAX_CHANGES_IN_REGION_RATIO * (read_length - kmer_length)),
                    modifications,
                    1,
                    checked_changes
                    );
        }
    }
        // second_kmer is not solid
    else {
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // not equal to the original character
            if (sequence_modified[kmer_length] != NEOCLEOTIDE[it_alter]) {
                // make a change
                second_kmer[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // new second_kmer is solid
                float kmer_quality;
                if (query_text(second_kmer, kmer_quality) == true) {
                    // if this k-mer is the last k-mer in a read
                    // running extend_a_kmer_3_prime_end is not needed any more
                    if ((read_length - kmer_length) == 1) {
                        // generate a new path
                        C_candidate_path new_candidate_path(candidate_path_in);

                        Single_mod pair_tmp(kmer_length, NEOCLEOTIDE[it_alter]);

                        new_candidate_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
                        new_candidate_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
                        new_candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

                        new_candidate_path.modified_bases.push_back(pair_tmp);
                        candidate_path_vector_tmp.push_back(new_candidate_path);
                    }
                    else if ((read_length - kmer_length) > 1) {
                        // trace  this kmer recursively and update candidate_path_vector_tmp

                        C_modification_with_quality modifications[max_read_length];
                        modifications[0].modification = NEOCLEOTIDE[it_alter];
                        modifications[0].quality = kmer_quality;

                        std::size_t checked_changes = 0;
                        extend_a_kmer_3_prime_end(
                                second_kmer,
                                1,
                                candidate_path_in,
                                candidate_path_vector_tmp,
                                (MAX_CHANGES_IN_REGION_RATIO * (read_length - kmer_length)) - 1,
                                modifications,
                                1,
                                checked_changes
                                );
                    }
                }
            }
        }
    }

    // check the solidness of the rightmost k-mers of each modified base
    std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

    // each candidate path
    for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); ++it_path) {
        std::string sequence_tmp(sequence_modified);
        perform_extend_out_right(sequence_tmp, *it_path, candidate_path_vector_tmp_tmp);
    }

    std::vector<C_candidate_path>::iterator best_it_path = choose_best_correction(candidate_path_vector_tmp_tmp);

    if (best_it_path != candidate_path_vector_tmp_tmp.end()) {
        correction_success = true;
        candidate_path_in = *best_it_path;
    }
}



//----------------------------------------------------------------------
// Modifies read until reaches a specified position.
//----------------------------------------------------------------------

inline void C_correct_read::extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector) {
    // generate a new k-mer
    std::string kmer_new(kmer.substr(1, kmer_length - 1));
    kmer_new.push_back(sequence_modified[index_kmer + kmer_length]);

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (query_text(kmer_new, kmer_quality) == true) {
        // if this k-mer is the last k-mer that can be modified
        // running extend_a_kmer_right is not needed any more
        current_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
        current_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
        current_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

        if ((index_kmer + 1) == index_last_mod) {
            candidate_path_vector.push_back(current_path);
        }
        else {
            extend_a_kmer(
                    kmer_new,
                    index_kmer + 1,
                    index_last_mod,
                    current_path,
                    candidate_path_vector
                    );
        }
    }
    else {
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // not equal to the original character
            if (sequence_modified[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
                // make a change
                kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query_text(kmer_new, kmer_quality) == true) {
                    // generate a new path
                    C_candidate_path temporary_path(current_path);
                    temporary_path.kmers_quality += kmer_quality;
#ifdef USE_KMER_MEDIAN
                    temporary_path.covering_kmers_weight_vector.push_back(kmer_quality);
#else
                    temporary_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif

                    Single_mod pair_tmp(index_kmer + kmer_length, NEOCLEOTIDE[it_alter]);

                    temporary_path.modified_bases.push_back(pair_tmp);

                    // if this k-mer is the last k-mer that can be modified
                    // running extend_a_kmer_right is not needed any more
                    if ((index_kmer + 1) == index_last_mod) {
                        candidate_path_vector.push_back(temporary_path);
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer(
                                kmer_new,
                                index_kmer + 1,
                                index_last_mod,
                                temporary_path,
                                candidate_path_vector
                                );
                    }
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Checks if the read can be extended towards 5'.
//----------------------------------------------------------------------

void C_correct_read::perform_extend_out_left(std::string& sequence_tmp, C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector_tmp_tmp) {
    // index_smallest_modified
    std::size_t index_smallest_modified(0);
    if (candidate_path.modified_bases.size() > 0) {
        index_smallest_modified = candidate_path.modified_bases[candidate_path.modified_bases.size() - 1].first;
    }

    // number of bases that should be extended
    std::size_t extend_amount;

    // calculate extend_amount
    // no extension is needed
    // kmer_length = 11, max_extension = 5
    // |0|0|0|0|0|0|0|0|0|0|1|1|1|-
    // |0|1|2|3|4|5|6|7|8|9|0|1|2|-
    // |<------------------->|      k = 11
    // |--------------------------- read
    //                     |<------ index_smallest_modified >= 10
    if (index_smallest_modified >= kmer_length - 1) {
        candidate_path_vector_tmp_tmp.push_back(candidate_path);
    }
        // extension is needed
    else {
        // applied the modified bases to sequence_tmp
        for (std::size_t it_base = 0; it_base < candidate_path.modified_bases.size(); it_base++) {
            // modify sequence_tmp
            sequence_tmp[candidate_path.modified_bases[it_base].first] = candidate_path.modified_bases[it_base].second;
        }

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

        // generate an initial k-mer
        std::string kmer_initial(sequence_tmp.substr(0, kmer_length - 1));
        kmer_initial = '0' + kmer_initial;

        float max_extended_kmer_quality = 0.0f;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[0] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            float kmer_quality;
            if (query_text(kmer_initial, kmer_quality) == true) {
                // if extend_amount == 1
                // running extend_out_left is not needed any more
                max_extended_kmer_quality = std::max(max_extended_kmer_quality, kmer_quality);
                if (extend_amount == 1) {
                    extension_success = true;
                }
                else if (!extension_success) {
                    // trace  this kmer recursively and update candidate_path_vector_tmp
                    extend_out_left(
                            kmer_initial,
                            1,
                            extend_amount,
                            extension_success
                            );
                }
            }
        }

        if (extension_success == true) {
            candidate_path.kmers_quality += max_extended_kmer_quality * EXTENSION_KMERS_WEIGHT;
#ifdef USE_KMER_MEDIAN
            it_path->covering_kmers_weight_vector.push_back(max_extended_kmer_quality * EXTENSION_KMERS_WEIGHT);
#else
            candidate_path.covering_kmers_weight += EXTENSION_KMERS_WEIGHT;
#endif
            candidate_path_vector_tmp_tmp.push_back(candidate_path);
        }
    }
}



//----------------------------------------------------------------------
// Checks if the read can be extended towards 3'.
//----------------------------------------------------------------------

void C_correct_read::perform_extend_out_right(std::string& sequence_tmp, C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector_tmp_tmp) {
    // index_largest_modified
    std::size_t index_largest_modified(read_length - 1);
    if (candidate_path.modified_bases.size() > 0) {
        index_largest_modified = candidate_path.modified_bases[candidate_path.modified_bases.size() - 1].first;
    }

    // number of bases that should be extended
    std::size_t extend_amount;

    // calculate extend_amount
    // no extension is needed
    // sequence.length() = 20, kmer_length = 11, max_extension = 5
    // |0|0|0|1|1|1|1|1|1|1|1|1|1|
    // |7|8|9|0|1|2|3|4|5|6|7|8|9|
    //     |<------------------->| k = 11
    // --------------------------| read
    // ----->|                     index_largest_modified <= 9
    if (index_largest_modified <= read_length - kmer_length) {
        candidate_path_vector_tmp_tmp.push_back(candidate_path);
    }
        // extension is needed
    else {
        // applied the modified bases to sequence_tmp
        for (std::size_t it_base = 0; it_base < candidate_path.modified_bases.size(); it_base++) {
            // modify sequence_tmp
            sequence_tmp[candidate_path.modified_bases[it_base].first] = (candidate_path).modified_bases[it_base].second;
        }

        // determine the number of extensions
        // sequence.length() = 20, kmer_length = 11, max_extension = 5
        // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
        // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
        //     |<------------------->|           k = 11
        // --------------------------|           read
        //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
        //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
        if (index_largest_modified <= read_length + max_extension - kmer_length) {
            extend_amount = kmer_length - (read_length - index_largest_modified);
        }
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //                 |<------->|           index_largest_modified > 15
        else {
            extend_amount = max_extension;
        }

        bool extension_success(false);

        // generate an initial k-mer
        // sequence.length() = 20, kmer_length = 11
        // |0|0|0|1|1|1|1|1|1|1|1|1|1|
        // |7|8|9|0|1|2|3|4|5|6|7|8|9|
        //       |<----------------->| kmer_length - 1 = 10
        // --------------------------| read
        //       |-|                   20 - 11 + 1 = 10
        std::string kmer_initial(sequence_tmp.substr(read_length - kmer_length + 1, kmer_length - 1));
        kmer_initial = kmer_initial + '0';

        float max_extended_kmer_quality = 0.0f;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            float kmer_quality;
            if (query_text(kmer_initial, kmer_quality) == true) {
                max_extended_kmer_quality = std::max(max_extended_kmer_quality, kmer_quality);
                // if extend_amount == 1
                // running extend_out_right is not needed any more
                if (extend_amount == 1) {
                    extension_success = true;
                }
                else if (!extension_success) {
                    // trace  this kmer recursively and update candidate_path_vector_tmp
                    extend_out_right(
                            kmer_initial,
                            1,
                            extend_amount,
                            extension_success
                            );
                }
            }
        }

        if (extension_success == true) {
            candidate_path.kmers_quality += max_extended_kmer_quality * EXTENSION_KMERS_WEIGHT;
#ifdef USE_KMER_MEDIAN
            candidate_path.covering_kmers_weight_vector.push_back(max_extended_kmer_quality * EXTENSION_KMERS_WEIGHT);
#else
            candidate_path.covering_kmers_weight += EXTENSION_KMERS_WEIGHT;
#endif
            candidate_path_vector_tmp_tmp.push_back(candidate_path);
        }
    }
}



//----------------------------------------------------------------------
// Proceeds correction towards 5' end.
//----------------------------------------------------------------------

inline void C_correct_read::extend_a_kmer_5_prime_end(const std::string& kmer, const std::size_t index_kmer, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t max_remaining_changes, C_modification_with_quality modifications[], std::size_t nesting, std::size_t& checked_changes) {
    if (candidate_path_vector.size() > MAX_EXTEND_CORRECTION_PATHS) {
        return;
    }
    // generate a new k-mer
    std::string kmer_new(kmer.substr(0, kmer_length - 1));
    kmer_new = sequence_modified[index_kmer - 1] + kmer_new;

    const bool is_low_quality_base = (quality_score[index_kmer - 1] - quality_score_offset < QS_CUTOFF);

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (!is_low_quality_base && query_text(kmer_new, kmer_quality) == true) {
        modifications[nesting].quality = kmer_quality;
        modifications[nesting].modification = kmer_new[0];

        // if this k-mer is the first k-mer in a read
        // running extend_a_kmer_5_prime_end is not needed any more
        if ((index_kmer - 1) == 0) {
            C_candidate_path candidate_path;
            for (std::size_t it = 0; it <= nesting; ++it) {
                if (modifications[it].modification != sequence_modified[nesting - it]) {
                    Single_mod single_modification(nesting - it, modifications[it].modification);
                    candidate_path.modified_bases.push_back(single_modification);
                }
                candidate_path.kmers_quality += modifications[it].quality;
#ifdef USE_KMER_MEDIAN
                candidate_path.covering_kmers_weight_vector.push_back(modifications[it].quality);
#else
                candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif
            }
            candidate_path_vector.push_back(candidate_path);
        }
        else if ((index_kmer - 1) > 0) {
            extend_a_kmer_5_prime_end(
                    kmer_new,
                    index_kmer - 1,
                    candidate_path_vector,
                    max_remaining_changes,
                    modifications,
                    nesting + 1,
                    checked_changes
                    );
        }
    }
    else {
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // not equal to the original character
            if (sequence_modified[index_kmer - 1] != NEOCLEOTIDE[it_alter] || is_low_quality_base) {
                // make a change
                kmer_new[0] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query_text(kmer_new, kmer_quality) == true) {
                    modifications[nesting].quality = kmer_quality;
                    modifications[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the first k-mer in a read
                    // running extend_a_kmer_5_prime_end is not needed any more
                    if ((index_kmer - 1) == 0) {
                        C_candidate_path candidate_path;
                        for (std::size_t it = 0; it <= nesting; ++it) {
                            if (modifications[it].modification != sequence_modified[nesting - it]) {
                                Single_mod single_modification(nesting - it, modifications[it].modification);
                                candidate_path.modified_bases.push_back(single_modification);
                            }
                            candidate_path.kmers_quality += modifications[it].quality;
#ifdef USE_KMER_MEDIAN
                            candidate_path.covering_kmers_weight_vector.push_back(modifications[it].quality);
#else
                            candidate_path.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif
                        }
                        candidate_path_vector.push_back(candidate_path);
                    }
                    else if ((index_kmer - 1) > 0
#ifdef LIMIT_MODIFICATIONS
                            && (max_remaining_changes > 1 || sequence_modified[index_kmer - 1] == NEOCLEOTIDE[it_alter])
#endif
                            ) {

                        ++checked_changes;
                        if (checked_changes > CHECK_MAX_CHANGES) {
                            return;
                        }

                        const std::size_t max_remaining_changes_new = (sequence_modified[index_kmer - 1] != NEOCLEOTIDE[it_alter] ? max_remaining_changes - 1 : max_remaining_changes);
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_5_prime_end(
                                kmer_new,
                                index_kmer - 1,
                                candidate_path_vector,
                                max_remaining_changes_new,
                                modifications,
                                nesting + 1,
                                checked_changes
                                );
                    }
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Proceeds correction towards 3' end.
//----------------------------------------------------------------------

inline void C_correct_read::extend_a_kmer_3_prime_end(const std::string& kmer, const std::size_t index_kmer, C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t max_remaining_changes, C_modification_with_quality modifications[], std::size_t nesting, std::size_t& checked_changes) {
    if (candidate_path_vector.size() > MAX_EXTEND_CORRECTION_PATHS) {
        return;
    }
    // generate a new k-mer
    std::string kmer_new(kmer.substr(1, kmer_length - 1));
    kmer_new = kmer_new + sequence_modified[index_kmer + kmer_length];

    const bool is_low_quality_base = (quality_score[index_kmer + kmer_length] - quality_score_offset < QS_CUTOFF);

    // kmer_new is a solid k-mer
    float kmer_quality;
    if (!is_low_quality_base && query_text(kmer_new, kmer_quality) == true) {
        modifications[nesting].quality = kmer_quality;
        modifications[nesting].modification = kmer_new[kmer_length - 1];

        // if this k-mer is the last k-mer in a read
        // running extend_a_kmer_3_prime_end is not needed any more
        if ((index_kmer + 1) == (read_length - kmer_length)) {
            C_candidate_path candidate_path_new(candidate_path);
            for (std::size_t it = 0; it <= nesting; ++it) {
                if (modifications[it].modification != sequence_modified[read_length - nesting + it - 1]) {
                    Single_mod single_modification(read_length - nesting + it - 1, modifications[it].modification);
                    candidate_path_new.modified_bases.push_back(single_modification);
                }
                candidate_path_new.kmers_quality += modifications[it].quality;
#ifdef USE_KMER_MEDIAN
                candidate_path_new.covering_kmers_weight_vector.push_back(modifications[it].quality);
#else
                candidate_path_new.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif
            }
            candidate_path_vector.push_back(candidate_path_new);
        }
        else if ((index_kmer + 1) < (read_length - kmer_length)) {
            extend_a_kmer_3_prime_end(
                    kmer_new,
                    index_kmer + 1,
                    candidate_path,
                    candidate_path_vector,
                    max_remaining_changes,
                    modifications,
                    nesting + 1,
                    checked_changes
                    );
        }
    }
    else {
        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // not equal to the original character
            if (sequence_modified[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter] || is_low_quality_base) {
                // make a change
                kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                // kmer_new is solid
                float kmer_quality;
                if (query_text(kmer_new, kmer_quality) == true) {
                    modifications[nesting].quality = kmer_quality;
                    modifications[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the last k-mer in a read
                    // running extend_a_kmer_3_prime_end is not needed any more
                    if ((index_kmer + 1) == (read_length - kmer_length)) {
                        C_candidate_path candidate_path_new(candidate_path);
                        for (std::size_t it = 0; it <= nesting; ++it) {
                            if (modifications[it].modification != sequence_modified[read_length - nesting + it - 1]) {
                                Single_mod single_modification(read_length - nesting + it - 1, modifications[it].modification);
                                candidate_path_new.modified_bases.push_back(single_modification);
                            }
                            candidate_path_new.kmers_quality += modifications[it].quality;
#ifdef USE_KMER_MEDIAN
                            candidate_path_new.covering_kmers_weight_vector.push_back(modifications[it].quality);
#else
                            candidate_path_new.covering_kmers_weight += COVERING_KMERS_WEIGHT;
#endif
                        }
                        candidate_path_vector.push_back(candidate_path_new);
                    }
                    else if ((index_kmer + 1) < (read_length - kmer_length)
#ifdef LIMIT_MODIFICATIONS
                            && (max_remaining_changes > 1 || sequence_modified[index_kmer - 1] == NEOCLEOTIDE[it_alter])
#endif
                            ) {
                        ++checked_changes;
                        if (checked_changes > CHECK_MAX_CHANGES) {
                            return;
                        }

                        // trace  this kmer recursively and update candidate_path_vector
                        const std::size_t max_remaining_changes_new = (sequence_modified[index_kmer - 1] != NEOCLEOTIDE[it_alter] ? max_remaining_changes - 1 : max_remaining_changes);
                        extend_a_kmer_3_prime_end(
                                kmer_new,
                                index_kmer + 1,
                                candidate_path,
                                candidate_path_vector,
                                max_remaining_changes_new,
                                modifications,
                                nesting + 1,
                                checked_changes
                                );
                    }
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Extends read towards 5'.
//----------------------------------------------------------------------

inline void C_correct_read::extend_out_left(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, bool& extension_success) {
    // generate a new k-mer
    std::string kmer_new(kmer.substr(0, kmer_length - 1));
    kmer_new = '0' + kmer_new;

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // generate kmer_new
        kmer_new[0] = NEOCLEOTIDE[it_alter];

        // kmer_new is solid
        float kmer_quality;
        if (query_text(kmer_new, kmer_quality) == true) {
            // if current num_extend = extend_amount
            // running extend_out_left is not needed any more
            if ((num_extend + 1) == extend_amount) {
                extension_success = true;
                break;
            }
            else {
                // trace  this kmer recursively
                extend_out_left(
                        kmer_new,
                        num_extend + 1,
                        extend_amount,
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
// Extends read towards 5'.
//----------------------------------------------------------------------

inline void C_correct_read::extend_out_right(const std::string& kmer, const std::size_t num_extend, const std::size_t extend_amount, bool& extension_success) {
    // generate a new k-mer
    std::string kmer_new(kmer.substr(1, kmer_length - 1));
    kmer_new = kmer_new + '0';

    // each alternative neocletide
    for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
        // generate kmer_new
        kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

        // kmer_new is solid
        float kmer_quality;
        if (query_text(kmer_new, kmer_quality) == true) {
            // if current num_extend = extend_amount
            // running extend_out_right is not needed any more
            if ((num_extend + 1) == extend_amount) {
                extension_success = true;
                break;
            }
            else {
                // trace  this kmer recursively
                extend_out_right(
                        kmer_new,
                        num_extend + 1,
                        extend_amount,
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

double C_correct_read::convert_quality_to_probability(char c) {
    c -= quality_score_offset;
    return std::pow(10.0, -c / 10.0);
}


//----------------------------------------------------------------------
// Checks in the KMC database if k-mer exists, if yes it returns kmer counter.
//----------------------------------------------------------------------

inline bool C_correct_read::query_text(const std::string& kmer, float& kmer_quality) {
    if (!kmer_api.from_string(kmer)) {
        return false;
    }

    if (kmc_file.CheckKmer(kmer_api, kmer_quality)) {
        return true;
    }

    kmer_api.reverse();

    return kmc_file.CheckKmer(kmer_api, kmer_quality);
}



//----------------------------------------------------------------------
// Chooses the best correction path.
//----------------------------------------------------------------------

std::vector<C_candidate_path>::iterator C_correct_read::choose_best_correction(std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path_vector.size() > 1) {
        // each path
        std::vector<C_candidate_path>::iterator it_path;
        std::vector<C_candidate_path>::iterator best_it_path;

        double best_kmer_quality = 0.0;

        // each candidate path
        for (it_path = candidate_path_vector.begin(); it_path != candidate_path_vector.end(); it_path++) {
            // each modification
            std::size_t it_first_mod = 0;
            double nucleotides_probability = 1.0;
            for (; it_first_mod < (*it_path).modified_bases.size(); it_first_mod++) {
                if (sequence[(*it_path).modified_bases[it_first_mod].first] != (*it_path).modified_bases[it_first_mod].second) {
                    nucleotides_probability = convert_quality_to_probability(quality_score[(*it_path).modified_bases[it_first_mod].first]);
                    break;
                }
            }
            for (std::size_t it_mod = it_first_mod + 1; it_mod < (*it_path).modified_bases.size(); it_mod++) {
                // add quality scores of modified bases
                if (sequence[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
                    nucleotides_probability *= convert_quality_to_probability(quality_score[(*it_path).modified_bases[it_mod].first]);
                }
            }

#ifdef USE_KMER_MEDIAN
            it_path->kmers_quality = median(it_path->covering_kmers_weight_vector);
#else
            it_path->kmers_quality /= it_path->covering_kmers_weight;
#endif
            it_path->kmers_quality *= nucleotides_probability;

            if (it_path->kmers_quality > best_kmer_quality) {
                best_kmer_quality = it_path->kmers_quality;
                best_it_path = it_path;
            }
        }

        // correction succeeds
        if (best_kmer_quality > MIN_BEST_KMER_QUALITY) {
            return best_it_path;
        }
        else {
            return candidate_path_vector.end();
        }
    }
        // only one path
        // correction succeeds
    else if (candidate_path_vector.size() == 1) {
        return candidate_path_vector.begin();
    }

    return candidate_path_vector.end();
}



//----------------------------------------------------------------------
// Applies changes into the read.
//----------------------------------------------------------------------

void C_correct_read::modify_errors(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors) {
    std::vector<C_candidate_path>::iterator best_it_path = choose_best_correction(candidate_path_vector);
    if (best_it_path != candidate_path_vector.end()) {
        // each modification
        for (std::vector< Single_mod >::iterator it_base = (*best_it_path).modified_bases.begin(); it_base != (*best_it_path).modified_bases.end(); ++it_base) {
            // update sequence_modification
            sequence_modification[(*it_base).first] = (*it_base).second;
            sequence_modified[(*it_base).first] = (*it_base).second;
            num_corrected_errors++;
        }
    }
}




//----------------------------------------------------------------------
// Applies changes into the read with respect of modification position.
//----------------------------------------------------------------------

void C_correct_read::modify_errors_first_kmer(std::vector<C_candidate_path>& candidate_path_vector, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2) {
    std::vector<C_candidate_path>::iterator best_it_path = choose_best_correction(candidate_path_vector);
    if (best_it_path != candidate_path_vector.end()) {
        for (std::size_t it_base = 0; it_base < (*best_it_path).modified_bases.size(); it_base++) {
            // filter out the bases that are equal to the original ones
            if (sequence[(*best_it_path).modified_bases[it_base].first] != (*best_it_path).modified_bases[it_base].second) {
                sequence_modification[(*best_it_path).modified_bases[it_base].first] = (*best_it_path).modified_bases[it_base].second;
                sequence_modified[(*best_it_path).modified_bases[it_base].first] = (*best_it_path).modified_bases[it_base].second;

                if (candidate_path_vector[0].modified_bases[it_base].first < kmer_length) {
                    num_corrected_errors1++;
                }
                else {
                    num_corrected_errors2++;
                }
            }
        }
    }
}



//----------------------------------------------------------------------
// Finds median from k-mers qualities.
//----------------------------------------------------------------------
#ifdef USE_KMER_MEDIAN

float C_correct_read::median(std::vector<float>& kmer_qualities) {
    std::sort(kmer_qualities.begin(), kmer_qualities.end());
    if (kmer_qualities.size() % 2 == 1) {
        return kmer_qualities[kmer_qualities.size() / 2];
    }
    return (kmer_qualities[kmer_qualities.size() / 2] + kmer_qualities[kmer_qualities.size() / 2 - 1]) * 0.5f;
}
#endif


//----------------------------------------------------------------------
// Removes modifications from a correction path.
//----------------------------------------------------------------------

void C_candidate_path::clear_path() {
    modified_bases.clear();
    kmers_quality = 0.0;
#ifndef USE_KMER_MEDIAN   
    covering_kmers_weight = 0.0;
#endif
}
