/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#include <thread>
#include "correct_errors.hpp"
#include "FileReader.h"
#include "Log.h"
#include <kmc_api/kmc_file.h>
#include <algorithm>
#include <mutex>
#include <vector>



//----------------------------------------------------------------------
// Replaces Ns in a read
//----------------------------------------------------------------------

std::size_t C_correction_utilities::replace_Ns(std::string& read) const {
    // change sequences to upper case
    transform(read.begin(), read.end(), read.begin(), ::toupper);

    // count and substitute Ns other characters
    std::size_t number_of_Ns = 0;
    for (std::size_t it = 0; it < read.length(); ++it) {
        bool is_nucleotide = false;
        for (std::size_t it_nucl = 0; it_nucl < NUM_NEOCLEOTIDE; ++it_nucl) {
            if (read[it] == NEOCLEOTIDE[it_nucl]) {
                is_nucleotide = true;
                break;
            }
        }
        if (!is_nucleotide) {
            ++number_of_Ns;
            read[it] = SUBST_CHAR;
        }
    }

    return number_of_Ns;
}



//----------------------------------------------------------------------
// Signal availability of a new file to merge.
//----------------------------------------------------------------------

void C_merge_temporary_files::signal_file_availability(size_t file_number) {
    std::unique_lock<std::mutex> lck(waiting_chunks_mutex);

    waiting_chunks.push(file_number);
    available_chunks.notify_one();
}



//----------------------------------------------------------------------
// Signal, that all files to merge have been created.
//----------------------------------------------------------------------

void C_merge_temporary_files::signal_finish() {
    std::unique_lock<std::mutex> lck(waiting_chunks_mutex);

    merging_finished = true;
    available_chunks.notify_one();
}



//----------------------------------------------------------------------
// Body of merging and removing of temporary files thread.
//----------------------------------------------------------------------

void C_merge_temporary_files::operator()() {
    std::ofstream output_file;

    size_t next_chunk = 0;

    while (true) {
        { // for unique_lock
            std::unique_lock<std::mutex> lck(waiting_chunks_mutex);

            available_chunks.wait(lck, [this, next_chunk] {return (!waiting_chunks.empty() && waiting_chunks.top() == next_chunk) || (merging_finished && waiting_chunks.empty()); });

            if (merging_finished && waiting_chunks.empty()) {
                break;
            }
            else {
                waiting_chunks.pop();
            }
        }

        // The first temporary file will be the beginning of the output file.
        if (next_chunk == 0) {
            output_file.open(readFileData.output_name, std::ios_base::binary | std::ios_base::app);

            if (!output_file.is_open()) {
                c_err << std::endl << "ERROR: Cannot create output file " << readFileData.output_name << "." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        // Append one file.
        else {
            std::string current_correction_info_file_name = readFileData.getCorrectionInfoName(next_chunk);

            std::ifstream input_file(current_correction_info_file_name, std::ios_base::binary);
            if (!input_file.is_open()) {
                c_err << std::endl << "ERROR: Cannot open temporary file " << current_correction_info_file_name << std::endl << std::endl;
            }

            output_file << input_file.rdbuf();

            if (!output_file) {
                c_err << std::endl << "ERROR: Cannot write to file " << readFileData.output_name << ". Probably disk is full." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }

            input_file.close();
            remove_error_correction_info_file(next_chunk);
        }
        ++next_chunk;
    }
}



//----------------------------------------------------------------------
// Creates correction threads.
//----------------------------------------------------------------------

void C_correct_errors::correct_errors_in_a_file(CKMCFile& kmc_file, CKMCFile* kmc_long_file) {
    if (!fastqReader.openFile(read_file_data.input_name, read_file_data.type)) {
        c_err << std::endl << "ERROR: Cannot open input file" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    std::thread reader_thread(std::ref(fastqReader));
    std::thread merger_thread(std::ref(merge_files));

    std::vector<std::thread> threads(n_threads);

    for (size_t it = 0; it < threads.size(); ++it) {
        threads[it] = std::thread(std::ref(*this), std::ref(kmc_file), kmc_long_file);
    }

    for (size_t it = 0; it < threads.size(); ++it) {
        threads[it].join();
    }

    merge_files.signal_finish();
    merger_thread.join();
    reader_thread.join();

    c_log << "     Number of corrected reads                  : " << global_num_corrected_reads << std::endl;
    c_log << "     Number of corrected substitutions          : " << global_num_corrected_substs << std::endl;
    if (correct_indels) {
        c_log << "     Number of corrected insertions             : " << global_num_corrected_ins << std::endl;
        c_log << "     Number of corrected deletions              : " << global_num_corrected_dels << std::endl;
    }

    if (verbose) {
        if (correct_indels) {
            c_log << "     Number of corrected both ins and dels      : " << global_num_corrected_pairs << std::endl;
            c_log << "     Number of regions with single deletions    : " << global_num_single_deletions << std::endl;
        }
        c_log << "     Number of corrected errors (step 1-1)      : " << global_num_corrected_errors_step1_1 << std::endl;
        c_log << "     Number of corrected errors (step 1-2)      : " << global_num_corrected_errors_step1_2 << std::endl;
        c_log << "     Number of corrected errors (step 1-3)      : " << global_num_corrected_errors_step1_3 << std::endl;
        c_log << "     Number of corrected errors (step 2-1)      : " << global_num_corrected_errors_step2_1 << std::endl;
        c_log << "     Number of corrected errors (step 2-2)      : " << global_num_corrected_errors_step2_2 << std::endl;

        c_log << "     Number of first k-mer corrections          : " << global_num_first_kmer_corrections << std::endl;
        c_log << "     Number of first k-mer successes            : " << global_num_first_kmer_successes << std::endl;
        if (kmc_long_file != nullptr) {
            c_log << "     Number of corrected reads checked with long k-mers:            " << global_num_checked_with_long_kmer << std::endl;
            c_log << "     Number of corrected reads filtered out with long k-mers:       " << global_num_filtered_with_long_kmer << std::endl;
            c_log << "     Number of corrected reads instantly accepted with long k-mers: " << global_num_instantly_accepted_with_long_kmer << std::endl;
        }
    }
}



//----------------------------------------------------------------------
// Correction thread body.
//----------------------------------------------------------------------

void C_correct_errors::operator()(CKMCFile& kmc_file, CKMCFile* kmc_long_file) {
    if (correct_indels) {
        C_correct_read<true> correct_read(quality_score_offset, kmer_length, long_kmer_length, extend, read_file_data.input_name, kmc_file, kmc_long_file, accept_filtered_with_long_kmers);

        // perform correction while some chunks are available
        while (true) {
            ReadsChunk readsChunk;
            if (!fastqReader.getChunk(readsChunk)) {
                break;
            }

            correct_errors_in_a_chunk(kmc_file, readsChunk.getChunkNo(), correct_read, readsChunk);

            merge_files.signal_file_availability(readsChunk.getChunkNo());
        }

        { // for lock_guard
            std::lock_guard<std::mutex> global_num_corrected_reads_guard(global_num_corrected_reads_mutex);
            global_num_corrected_reads += correct_read.num_corrected_reads;

            global_num_corrected_errors_step1_1 += correct_read.num_corrected_errors_step1_1;
            global_num_corrected_errors_step1_2 += correct_read.num_corrected_errors_step1_2;
            global_num_corrected_errors_step1_3 += correct_read.num_corrected_errors_step1_3;
            global_num_corrected_errors_step2_1 += correct_read.num_corrected_errors_step2_1;
            global_num_corrected_errors_step2_2 += correct_read.num_corrected_errors_step2_2;

            global_num_corrected_substs += correct_read.num_corrected_substs;
            global_num_corrected_ins += correct_read.num_corrected_ins;
            global_num_corrected_dels += correct_read.num_corrected_dels;
            global_num_corrected_pairs += correct_read.num_corrected_pairs;
            global_num_first_kmer_corrections += correct_read.num_first_kmer_corrections;
            global_num_first_kmer_successes += correct_read.num_first_kmer_successes;
            global_num_single_deletions += correct_read.num_single_deletions;
            global_num_checked_with_long_kmer += correct_read.num_checked_with_long_kmer;
            global_num_filtered_with_long_kmer += correct_read.num_filtered_with_long_kmer;
            global_num_instantly_accepted_with_long_kmer += correct_read.num_instantly_accepted_with_long_kmer;
            global_num_oriented_reads += correct_read.num_oriented_reads;
        }
    }
    else {
        C_correct_read<false> correct_read(quality_score_offset, kmer_length, long_kmer_length, extend, read_file_data.input_name, kmc_file, kmc_long_file, accept_filtered_with_long_kmers);

        // perform correction while some chunks are available
        while (true) {
            ReadsChunk readsChunk;
            if (!fastqReader.getChunk(readsChunk)) {
                break;
            }

            correct_errors_in_a_chunk(kmc_file, readsChunk.getChunkNo(), correct_read, readsChunk);

            merge_files.signal_file_availability(readsChunk.getChunkNo());
        }

        { // for lock_guard
            std::lock_guard<std::mutex> global_num_corrected_reads_guard(global_num_corrected_reads_mutex);
            global_num_corrected_reads += correct_read.num_corrected_reads;

            global_num_corrected_errors_step1_1 += correct_read.num_corrected_errors_step1_1;
            global_num_corrected_errors_step1_2 += correct_read.num_corrected_errors_step1_2;
            global_num_corrected_errors_step1_3 += correct_read.num_corrected_errors_step1_3;
            global_num_corrected_errors_step2_1 += correct_read.num_corrected_errors_step2_1;
            global_num_corrected_errors_step2_2 += correct_read.num_corrected_errors_step2_2;

            global_num_corrected_substs += correct_read.num_corrected_substs;
            global_num_first_kmer_corrections += correct_read.num_first_kmer_corrections;
            global_num_first_kmer_successes += correct_read.num_first_kmer_successes;
            global_num_checked_with_long_kmer += correct_read.num_checked_with_long_kmer;
            global_num_filtered_with_long_kmer += correct_read.num_filtered_with_long_kmer;
            global_num_instantly_accepted_with_long_kmer += correct_read.num_instantly_accepted_with_long_kmer;
            global_num_oriented_reads += correct_read.num_oriented_reads;
        }
    }
}




//----------------------------------------------------------------------
// Saves corrected read into a temporary file.
//----------------------------------------------------------------------

void C_correct_errors::write_a_read_to_fastq(const std::string& header, const std::string& connector, const std::string& qualities, const std::string& sequence, FileReader& f_error_correction) {
    f_error_correction.putString(header);
    f_error_correction.putString("\n");

    f_error_correction.putString(sequence);
    f_error_correction.putString("\n");

    f_error_correction.putString(connector);
    f_error_correction.putString("\n");

    f_error_correction.putString(qualities);
    f_error_correction.putString("\n");
}



//----------------------------------------------------------------------
// Determines if a header or connector contains information about read length. Potentially changes it.
//----------------------------------------------------------------------

template<>
void C_correct_errors::change_headers_length_tag<false>(std::string& header, std::size_t original_length, std::size_t new_length) {
    // do nothing
}



template<>
void C_correct_errors::change_headers_length_tag<true>(std::string& header, std::size_t original_length, std::size_t new_length) {
    if (original_length == new_length) {
        return;
    }

    static const std::string tag_seq = " length=";
    const std::size_t tag_pos = header.rfind(tag_seq);
    if (tag_pos == std::string::npos) {
        return;
    }

    const std::string original_length_seq = std::to_string(original_length);

    // check, if the length is originally properly given
    if (header.substr(tag_pos + tag_seq.length()) == original_length_seq) {
        const std::string new_length_seq = std::to_string(new_length);
        header.replace(tag_pos + tag_seq.length(), original_length_seq.length(), new_length_seq);
    }
}
