/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.1
 * 
 */

#include "correct_read.hpp"



 //----------------------------------------------------------------------
 // Redetermines path rate after extending the read.
 //----------------------------------------------------------------------

void C_candidate_path::attach_extending_rate(double max_extending_kmer_quality) {
    saved_kmers_quality += max_extending_kmer_quality * EXTENSION_KMERS_WEIGHT;
    saved_covering_kmers_weight += EXTENSION_KMERS_WEIGHT;
    rate_internal(saved_kmers_quality, saved_covering_kmers_weight);
}



 //----------------------------------------------------------------------
 // Attaches the given correction path at the beginning.
 //----------------------------------------------------------------------

void C_candidate_path::attach_modifications_at_beginning(const C_candidate_path& beginning_path) {
    std::vector<Single_mod> new_modified_bases(beginning_path.modifications);
    new_modified_bases.reserve(new_modified_bases.size() + modifications.size());
    new_modified_bases.insert(new_modified_bases.end(), modifications.begin(), modifications.end());
    modifications.swap(new_modified_bases);

    num_insertions += beginning_path.num_insertions;
    read_length_change += beginning_path.read_length_change;
    beginning_attached_path_modifications = static_cast<int>(beginning_path.modifications.size());

    for (std::size_t it = beginning_attached_path_modifications; it < modifications.size(); ++it) {
        modifications[it].pos += beginning_path.read_length_change;
    }

    saved_kmers_quality += beginning_path.saved_kmers_quality;
    saved_covering_kmers_weight += beginning_path.saved_covering_kmers_weight;
    nucleotides_probability *= beginning_path.nucleotides_probability;

    rate_internal(saved_kmers_quality, saved_covering_kmers_weight);
}



//----------------------------------------------------------------------
// Determines the path rate.
//----------------------------------------------------------------------

void C_candidate_path::rate_internal(double kmers_quality, double covering_kmers_weight) {
    if (num_insertions > 0 && num_insertions == modifications.size()) {
        path_rate = INSERTION_RATE;
    }
    else {
        path_rate = kmers_quality;
        path_rate /= covering_kmers_weight;
        path_rate *= nucleotides_probability;
    }
    saved_kmers_quality = kmers_quality;
    saved_covering_kmers_weight = covering_kmers_weight;
}



void C_candidate_path::rate(double kmers_quality, std::size_t num_covering_kmers)
{
    rate_internal(kmers_quality * COVERING_KMERS_WEIGHT, num_covering_kmers * COVERING_KMERS_WEIGHT);
}



//----------------------------------------------------------------------
// Modifies read until reaches a specified position.
//----------------------------------------------------------------------

template<>
void C_correct_read<true>::extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid) {
    assert(index_kmer <= index_last_mod);
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
    double kmer_quality;
    if (!is_low_quality_base && query.query_text(kmer_new, kmer_quality) == true) {
        modifications_sequence_stack[nesting].quality = kmer_quality;
        modifications_sequence_stack[nesting].modification = kmer_new.back();

        // if this k-mer is the last k-mer that can be modified
        // running extend_a_kmer_right is not needed any more
        if (index_kmer == index_last_mod) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                return;
            }
            create_modification_path_with_indels_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting);
        }
        else {
            extend_a_kmer(
                kmer_new,
                index_kmer + 1,
                index_last_mod,
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
                double kmer_quality;
                if (query.query_text(kmer_new, kmer_quality) == true) {
                    solid_found = true;

                    modifications_sequence_stack[nesting].quality = kmer_quality;
                    modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the last k-mer that can be modified
                    // running extend_a_kmer_right is not needed any more
                    if (index_kmer == index_last_mod) {
                        if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                            return;
                        }
                        if (!create_modification_path_with_indels_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting)) {
                            if (!correct_last_deletion(kmer_new, candidate_path_vector, index_kmer, nesting + 1, checked_changes)) {
                                check_is_new_potential_indel(nesting, index_kmer + kmer_length - 1, true);
                            }
                        }
                    }
                    else {
                        // trace this kmer recursively and update candidate_path_vector
                        extend_a_kmer(
                            kmer_new,
                            index_kmer + 1,
                            index_last_mod,
                            candidate_path_vector,
                            nesting + 1,
                            checked_changes,
                            max_remaining_changes - 1,
                            max_remaining_non_solid
                            );

                        check_is_new_potential_indel(nesting, index_kmer + kmer_length - 1);
                    }
                }
            }
        }

        // try to accept non-solid k-mer
        if (!solid_found) {
            if (max_remaining_non_solid > 0) {
                kmer_new.back() = sequence_modified[index_kmer + kmer_length - 1];

                modifications_sequence_stack[nesting].quality = 0;
                modifications_sequence_stack[nesting].modification = kmer_new.back();

                // if this k-mer is the last k-mer in a read
                // running extend_a_kmer_3_prime_end is not needed any more
                if (index_kmer == index_last_mod) {
                    if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                        return;
                    }
                    if (!create_modification_path_with_indels_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting)) {
                        check_is_new_potential_indel(nesting, index_kmer + kmer_length - 1, true);
                    }
                }
                else {
                    extend_a_kmer(
                        kmer_new,
                        index_kmer + 1,
                        index_last_mod,
                        candidate_path_vector,
                        nesting + 1,
                        checked_changes,
                        max_remaining_changes,
                        max_remaining_non_solid - 1
                        );
                }
            }
        }
        correct_indel(kmer_new, index_kmer, index_last_mod, candidate_path_vector, nesting, checked_changes, max_remaining_changes, max_remaining_non_solid);
    }
}



template<>
void C_correct_read<false>::extend_a_kmer(const std::string& kmer, const std::size_t index_kmer, const std::size_t index_last_mod, std::vector<C_candidate_path>& candidate_path_vector, std::size_t nesting, std::size_t& checked_changes, const int max_remaining_changes, const std::size_t max_remaining_non_solid) {
    assert(index_kmer <= index_last_mod);
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
    double kmer_quality;
    if (!is_low_quality_base && query.query_text(kmer_new, kmer_quality) == true) {
        modifications_sequence_stack[nesting].quality = kmer_quality;
        modifications_sequence_stack[nesting].modification = kmer_new.back();

        // if this k-mer is the last k-mer that can be modified
        // running extend_a_kmer_right is not needed any more
        if (index_kmer == index_last_mod) {
            if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                return;
            }
            create_modification_path_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting);
        }
        else {
            extend_a_kmer(
                kmer_new,
                index_kmer + 1,
                index_last_mod,
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
                double kmer_quality;
                if (query.query_text(kmer_new, kmer_quality) == true) {
                    solid_found = true;

                    modifications_sequence_stack[nesting].quality = kmer_quality;
                    modifications_sequence_stack[nesting].modification = NEOCLEOTIDE[it_alter];

                    // if this k-mer is the last k-mer that can be modified
                    // running extend_a_kmer_right is not needed any more
                    if (index_kmer == index_last_mod) {
                        if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                            return;
                        }
                        create_modification_path_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting);
                    }
                    else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer(
                            kmer_new,
                            index_kmer + 1,
                            index_last_mod,
                            candidate_path_vector,
                            nesting + 1,
                            checked_changes,
                            max_remaining_changes - 1,
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
                if (index_kmer == index_last_mod) {
                    if (candidate_path_vector.size() >= MAX_EXTEND_CORRECTION_PATHS) {
                        return;
                    }
                    create_modification_path_towards_3_prime_internal(candidate_path_vector, index_kmer, nesting);
                }
                else {
                    extend_a_kmer(
                        kmer_new,
                        index_kmer + 1,
                        index_last_mod,
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
// Checks if the read can be extended towards 5'.
//----------------------------------------------------------------------

template<>
bool C_correct_read<true>::perform_extend_out_left(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path.modifications.empty()) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // index_smallest_modified
    std::size_t index_smallest_modified(candidate_path.modifications.back().pos);

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
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // extension is needed
    else {
        std::string sequence_tmp(sequence_modified);
        // apply the modified bases to sequence_tmp
        apply_path_to_temporary_read(candidate_path, sequence_tmp);

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
        std::string kmer_initial(sequence_tmp.substr(0, kmer_length - 1));
        kmer_initial = '0' + kmer_initial;

        double max_extended_kmer_quality = 0.0f;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[0] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            double kmer_quality = 0.0f;
            bool solid_kmer = query.query_text(kmer_initial, kmer_quality);
            if (solid_kmer || max_remaining_non_solid > 0) {
                // if extend_amount == 1
                // running extend_out_left is not needed any more
                max_extended_kmer_quality = std::max(max_extended_kmer_quality, kmer_quality);
                if (extend_amount == 1) {
                    if (solid_kmer) {
                        extension_success = true;
                    }
                }
                else if (!extension_success) {
                    // trace this kmer recursively and update candidate_path_vector_tmp
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
            candidate_path.attach_extending_rate(max_extended_kmer_quality);

            candidate_path_vector.push_back(candidate_path);
        }
        return extension_success;
    }
}



template<>
bool C_correct_read<false>::perform_extend_out_left(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path.modifications.empty()) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // index_smallest_modified
    std::size_t index_smallest_modified(candidate_path.modifications.back().pos);

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
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // extension is needed
    else {
        std::string sequence_tmp(sequence_modified);
        // apply the modified bases to sequence_tmp
        apply_path_to_temporary_read(candidate_path, sequence_tmp);

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
        std::string kmer_initial(sequence_tmp.substr(0, kmer_length - 1));
        kmer_initial = '0' + kmer_initial;

        double max_extended_kmer_quality = 0.0f;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[0] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            double kmer_quality = 0.0f;
            bool solid_kmer = query.query_text(kmer_initial, kmer_quality);
            if (solid_kmer || max_remaining_non_solid > 0) {
                // if extend_amount == 1
                // running extend_out_left is not needed any more
                max_extended_kmer_quality = std::max(max_extended_kmer_quality, kmer_quality);
                if (extend_amount == 1) {
                    if (solid_kmer) {
                        extension_success = true;
                    }
                }
                else if (!extension_success) {
                    // trace this kmer recursively and update candidate_path_vector_tmp
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
            candidate_path.attach_extending_rate(max_extended_kmer_quality);

            candidate_path_vector.push_back(candidate_path);
        }
        return extension_success;
    }
}



//----------------------------------------------------------------------
// Checks if the read can be extended towards 3'.
//----------------------------------------------------------------------

template<>
bool C_correct_read<true>::perform_extend_out_right(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path.modifications.empty()) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }

    // index_largest_modified
    int index_largest_modified(candidate_path.modifications.back().pos);
    index_largest_modified += candidate_path.get_read_length_change();

    int last_kmer_pos = static_cast<int>(read_length - kmer_length + candidate_path.get_read_length_change());

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
    if (index_largest_modified <= last_kmer_pos) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // extension is needed
    else {
        std::string sequence_tmp(sequence_modified);
        // apply the modified bases to sequence_tmp
        apply_path_to_temporary_read(candidate_path, sequence_tmp);

        // determine the number of extensions
        // sequence.length() = 20, kmer_length = 11, max_extension = 5
        // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
        // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
        //     |<------------------->|           k = 11
        // --------------------------|           read
        //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
        //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
        if (index_largest_modified <= last_kmer_pos + static_cast<int>(max_extension)) {
            extend_amount = kmer_length - (last_kmer_pos + kmer_length - index_largest_modified);
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
        std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil(extend_amount * MAX_NON_SOLID));

        // generate an initial k-mer
        // sequence.length() = 20, kmer_length = 11
        // |0|0|0|1|1|1|1|1|1|1|1|1|1|
        // |7|8|9|0|1|2|3|4|5|6|7|8|9|
        //       |<----------------->| kmer_length - 1 = 10
        // --------------------------| read
        //       |-|                   20 - 11 + 1 = 10
        std::string kmer_initial(sequence_tmp.substr(last_kmer_pos + 1, kmer_length - 1));
        kmer_initial += '0';

        double max_extended_kmer_quality = 0.0f;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            double kmer_quality = 0.0f;
            bool solid_kmer = query.query_text(kmer_initial, kmer_quality);
            if (solid_kmer || max_remaining_non_solid > 0) {
                // if extend_amount == 1
                // running extend_out_right is not needed any more
                max_extended_kmer_quality = std::max(max_extended_kmer_quality, kmer_quality);
                if (extend_amount == 1) {
                    if (solid_kmer) {
                        extension_success = true;
                    }
                }
                else if (!extension_success) {
                    // trace this kmer recursively and update candidate_path_vector_tmp
                    extend_out_right(
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
            candidate_path.attach_extending_rate(max_extended_kmer_quality);

            candidate_path_vector.push_back(candidate_path);
        }
        return extension_success;
    }
}



template<>
bool C_correct_read<false>::perform_extend_out_right(C_candidate_path& candidate_path, std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path.modifications.empty()) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // index_largest_modified
    std::size_t index_largest_modified(candidate_path.modifications.back().pos);

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
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    // extension is needed
    else {
        std::string sequence_tmp(sequence_modified);
        // apply the modified bases to sequence_tmp
        apply_path_to_temporary_read(candidate_path, sequence_tmp);

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
        std::size_t max_remaining_non_solid = static_cast<std::size_t>(std::ceil(extend_amount * MAX_NON_SOLID));

        // generate an initial k-mer
        // sequence.length() = 20, kmer_length = 11
        // |0|0|0|1|1|1|1|1|1|1|1|1|1|
        // |7|8|9|0|1|2|3|4|5|6|7|8|9|
        //       |<----------------->| kmer_length - 1 = 10
        // --------------------------| read
        //       |-|                   20 - 11 + 1 = 10
        std::string kmer_initial(sequence_tmp.substr(read_length - kmer_length + 1, kmer_length - 1));
        kmer_initial += '0';

        double max_extended_kmer_quality = 0.0f;

        // each alternative neocletide
        for (std::size_t it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            double kmer_quality = 0.0f;
            bool solid_kmer = query.query_text(kmer_initial, kmer_quality);
            if (solid_kmer || max_remaining_non_solid > 0) {
                // if extend_amount == 1
                // running extend_out_right is not needed any more
                max_extended_kmer_quality = std::max(max_extended_kmer_quality, kmer_quality);
                if (extend_amount == 1) {
                    if (solid_kmer) {
                        extension_success = true;
                    }
                }
                else if (!extension_success) {
                    // trace this kmer recursively and update candidate_path_vector_tmp
                    extend_out_right(
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
            candidate_path.attach_extending_rate(max_extended_kmer_quality);

            candidate_path_vector.push_back(candidate_path);
        }
        return extension_success;
    }
}

//----------------------------------------------------------------------
// Extends corrected region towards 3'.
//----------------------------------------------------------------------

template<>
bool C_correct_read<true>::perform_extend_right_inner_region(C_candidate_path & candidate_path, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_last_mod_kmer) {
    if (candidate_path.modifications.empty()) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }

    index_last_mod_kmer += candidate_path.get_read_length_change();

    std::size_t index_last_kmer_to_check = candidate_path.modifications.back().pos - kmer_length + 1 + max_extension;

    if (index_last_kmer_to_check > index_last_mod_kmer) {
        // generate a temporary sequence
        std::string sequence_tmp(sequence_modified);
        apply_path_to_temporary_read(candidate_path, sequence_tmp);

        double max_extended_kmer_quality = 0.0f;

        // check k-mers
        std::size_t num_success(0);
        {
            // check the first extending k-mer to obtain its quality
            const char* current_kmer = sequence_tmp.c_str() + index_last_mod_kmer + 1; // add 1 to omit the last already checked k-mer

            if (query.query_text(current_kmer, max_extended_kmer_quality) == true) {
                num_success++;
            }
        }
        for (std::size_t it_check = index_last_mod_kmer + 2; it_check <= index_last_kmer_to_check; it_check++) {
            const char* current_kmer = sequence_tmp.c_str() + it_check;

            double kmer_quality;
            if (query.query_text(current_kmer, kmer_quality) == true) {
                num_success++;
            }
        }

        std::size_t max_non_solid = static_cast<std::size_t>(std::ceil((index_last_kmer_to_check - index_last_mod_kmer + 1) * MAX_NON_SOLID));

        if (num_success >= (index_last_kmer_to_check - index_last_mod_kmer + 1) - max_non_solid) {
            candidate_path.attach_extending_rate(max_extended_kmer_quality);

            candidate_path_vector.push_back(candidate_path);
            return true;
        }
    }
    // checking is not needed
    else {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    return false;
}



template<>
bool C_correct_read<false>::perform_extend_right_inner_region(C_candidate_path & candidate_path, std::vector<C_candidate_path>& candidate_path_vector, std::size_t index_last_mod_kmer) {
    if (candidate_path.modifications.empty()) {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }

    std::size_t index_last_kmer_to_check(candidate_path.modifications.back().pos - kmer_length + 1 + max_extension);

    if (index_last_kmer_to_check > index_last_mod_kmer) {
        // generate a temporary sequence
        std::string sequence_tmp(sequence_modified);
        apply_path_to_temporary_read(candidate_path, sequence_tmp);

        double max_extended_kmer_quality = 0.0f;

        // check k-mers
        std::size_t num_success(0);
        {
            // check the first extending k-mer to obtain its quality
            const char* current_kmer = sequence_tmp.c_str() + index_last_mod_kmer + 1; // add 1 to omit the last already checked k-mer

            if (query.query_text(current_kmer, max_extended_kmer_quality) == true) {
                num_success++;
            }
        }
        for (std::size_t it_check = index_last_mod_kmer + 2; it_check <= index_last_kmer_to_check; it_check++) { // add 1 to omit the last already checked k-mer
            const char* current_kmer = sequence_tmp.c_str() + it_check;

            double kmer_quality;
            if (query.query_text(current_kmer, kmer_quality) == true) {
                num_success++;
            }
        }

        std::size_t max_non_solid = static_cast<std::size_t>(std::ceil((index_last_kmer_to_check - index_last_mod_kmer + 1) * MAX_NON_SOLID));

        if (num_success >= (index_last_kmer_to_check - index_last_mod_kmer + 1) - max_non_solid) {
            candidate_path.attach_extending_rate(max_extended_kmer_quality);

            candidate_path_vector.push_back(candidate_path);
            return true;
        }
    }
    // checking is not needed
    else {
        candidate_path_vector.push_back(candidate_path);
        return true;
    }
    return false;
}



//----------------------------------------------------------------------
// Applies changes present in the path to the output sequence.
//----------------------------------------------------------------------

template<>
void C_correct_read<true>::apply_path_to_temporary_read(const C_candidate_path & candidate_path_in, std::string & sequence_out) {
    for (const Single_mod& single_mod : candidate_path_in.modifications) {
        assert(single_mod.pos < sequence_out.length() || (single_mod.pos == sequence_out.length() && single_mod.insertion)); // the condition after || is added just for completeness, actually it won't happen
        if (single_mod.insertion) {
            sequence_out.erase(single_mod.pos, 1);
        }
        else if (single_mod.deletion) {
            sequence_out.insert(single_mod.pos, 1, single_mod.modification);
        }
        else { // substitution
            sequence_out[single_mod.pos] = single_mod.modification;
        }
    }
}



template<>
void C_correct_read<false>::apply_path_to_temporary_read(const C_candidate_path & candidate_path_in, std::string & sequence_out) {
    for (const Single_mod& single_mod : candidate_path_in.modifications) {
        assert(single_mod.pos < sequence_out.length());
        // substitution
        sequence_out[single_mod.pos] = single_mod.modification;
    }
}



//----------------------------------------------------------------------
// Applies changes present in the path to the output sequence and qualities sequence.
//----------------------------------------------------------------------

template<>
void C_correct_read<true>::apply_path_to_read(const C_candidate_path & candidate_path_in, std::string & sequence_out, std::string& qualities_out, std::size_t& num_corrected_errors) {
    // optimization for oft-used situation, where path contains single substitution only
    if (candidate_path_in.modifications.size() == 1 && !candidate_path_in.modifications.front().insertion && !candidate_path_in.modifications.front().deletion) {
        sequence_out[candidate_path_in.modifications.front().pos] = candidate_path_in.modifications.front().modification;
        num_corrected_errors++;
    }
    else {
        bool was_ins = false, was_del = false;

        for (const Single_mod& single_mod : candidate_path_in.modifications) {
            if (single_mod.insertion) {
                sequence_out.erase(single_mod.pos, 1);
                qualities_out.erase(single_mod.pos, 1);
                was_ins = true;
                num_corrected_ins++;
            }
            else if (single_mod.deletion) {
                sequence_out.insert(single_mod.pos, 1, single_mod.modification);
                qualities_out.insert(single_mod.pos, 1, DUMMY_QUALITY_VALUE + quality_score_offset);
                was_del = true;
                num_corrected_dels++;
            }
            else { // substitution
                sequence_out[single_mod.pos] = single_mod.modification;
                num_corrected_substs++;
            }
            num_corrected_errors++;
        }

        if (was_ins && was_del) {
            num_corrected_pairs++;
        }
    }
}



template<>
void C_correct_read<false>::apply_path_to_read(const C_candidate_path & candidate_path_in, std::string & sequence_out, std::string& qualities_out, std::size_t& num_corrected_errors) {
    for (const Single_mod& single_mod : candidate_path_in.modifications) {
        // substitution
        sequence_out[single_mod.pos] = single_mod.modification;
        num_corrected_errors++;
        num_corrected_substs++;
    }
    // change of qualities_out is not needed
}



template<>
void C_correct_read<true>::apply_path_to_read(const C_candidate_path& candidate_path_in, std::string& sequence_out, std::string& qualities_out, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2) {
    bool was_ins = false, was_del = false;

    std::size_t mod_pos = 0;
    for (const Single_mod& single_mod : candidate_path_in.modifications) {
        if (single_mod.insertion) {
            sequence_out.erase(single_mod.pos, 1);
            qualities_out.erase(single_mod.pos, 1);
            was_ins = true;
            num_corrected_ins++;
        }
        else if (single_mod.deletion) {
            sequence_out.insert(single_mod.pos, 1, single_mod.modification);
            qualities_out.insert(single_mod.pos, 1, DUMMY_QUALITY_VALUE + quality_score_offset);
            was_del = true;
            num_corrected_dels++;
        }
        else { // substitution
            sequence_out[single_mod.pos] = single_mod.modification;
            num_corrected_substs++;
        }
        if (mod_pos < candidate_path_in.get_beginning_attached_path_modifications()) {
            num_corrected_errors1++;
        }
        else {
            num_corrected_errors2++;
        }
        mod_pos++;
    }

    if (was_ins && was_del) {
        num_corrected_pairs++;
    }
}



template<>
void C_correct_read<false>::apply_path_to_read(const C_candidate_path& candidate_path_in, std::string& sequence_out, std::string& qualities_out, std::size_t& num_corrected_errors1, std::size_t& num_corrected_errors2) {
    std::size_t mod_pos = 0;
    for (const Single_mod& single_mod : candidate_path_in.modifications) {
        // substitution
        sequence_out[single_mod.pos] = single_mod.modification;
        if (mod_pos < candidate_path_in.get_beginning_attached_path_modifications()) {
            num_corrected_errors1++;
        }
        else {
            num_corrected_errors2++;
        }
        num_corrected_substs++;
        mod_pos++;
    }
    // change of qualities_out is not needed
}



//----------------------------------------------------------------------
// Chooses the best correction path.
//----------------------------------------------------------------------

template<>
std::vector<C_candidate_path>::iterator C_correct_read<true>::choose_best_correction(std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path_vector.size() > 1) {
        // each path
        std::vector<C_candidate_path>::iterator it_path;
        std::vector<C_candidate_path>::iterator best_it_path;

        double best_rate = 0.0;

        // each candidate path
        for (it_path = candidate_path_vector.begin(); it_path != candidate_path_vector.end(); ++it_path) {
            double rate = it_path->get_path_rate();

            if (rate > best_rate) {
                best_rate = rate;
                best_it_path = it_path;
            }
        }

        // correction succeeds
        if (best_rate > MIN_RATE) {
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



template<>
std::vector<C_candidate_path>::iterator C_correct_read<false>::choose_best_correction(std::vector<C_candidate_path>& candidate_path_vector) {
    if (candidate_path_vector.size() > 1) {
        // each path
        std::vector<C_candidate_path>::iterator it_path;
        std::vector<C_candidate_path>::iterator best_it_path;

        double best_rate = 0.0;

        // each candidate path
        for (it_path = candidate_path_vector.begin(); it_path != candidate_path_vector.end(); ++it_path) {
            double rate = it_path->get_path_rate();

            if (rate > best_rate) {
                best_rate = rate;
                best_it_path = it_path;
            }
        }

        // correction succeeds
        if (best_rate > MIN_RATE) {
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
// Creates modification of a single k-mer vector.
//----------------------------------------------------------------------

template<>
bool C_correct_read<true>::create_modification_path_with_single_indel(std::vector<C_candidate_path>& candidate_path_vector, double kmer_quality, std::size_t pos, bool insertion /*= true*/, char modification /*= '0'*/) {
    C_candidate_path candidate_path;

    if (insertion) {
        candidate_path.add_insertion(pos);
    }
    else {
        candidate_path.add_deletion(pos, modification);
    }
    candidate_path.rate(kmer_quality, 1U);

    return perform_extend_out_first_kmer(candidate_path, candidate_path_vector);
}



//----------------------------------------------------------------------
// To allow calling of C_correct_read<true>::create_modification_path_with_single_indel from first k-mer correction, with no specialization of that function.
//----------------------------------------------------------------------

template<>
bool C_correct_read<false>::create_modification_path_with_single_indel(std::vector<C_candidate_path>& candidate_path_vector, double kmer_quality, std::size_t pos, bool insertion /*= true*/, char modification /*= '0'*/) {
    assert(false);
    return false;
}
