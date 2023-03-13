/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#ifndef _CORRECT_H
#define _CORRECT_H

#include "FastqReader.h"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "correct_read.hpp"
#include "Log.h"
#include <string>
#include <mutex>
#include <queue>
#include <condition_variable>



class C_correction_utilities {
public:
    std::size_t replace_Ns(std::string& read) const;
};



class C_merge_temporary_files {
private:
    C_log c_log;
    C_log c_err;

    C_correction_utilities correction_utilities;

    const ReadFileData readFileData;

    std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> waiting_chunks;
    std::mutex waiting_chunks_mutex;
    std::condition_variable available_chunks;

    bool merging_finished;

    void remove_error_correction_info_file(const std::size_t file_index) const {
        remove(readFileData.getCorrectionInfoName(file_index).c_str());
    }

public:
    C_merge_temporary_files(const ReadFileData& _readFileData) :
        c_log(std::cout),
        c_err(std::cerr),
        readFileData(_readFileData),
        merging_finished(false)
    {}

    void signal_file_availability(size_t file_number);
    void signal_finish();

    void operator()();
};




class C_correct_errors {
private:
    C_log c_log;
    C_log c_err;

    C_correction_utilities correctionUtilities;

    std::size_t quality_score_offset;
    std::size_t kmer_length;
    std::size_t long_kmer_length;
    std::size_t extend;
    bool accept_filtered_with_long_kmers;

    bool correct_indels;
    bool b_change_headers_length_tag;

    bool mark_corrected;

    ReadFileData read_file_data;


    C_merge_temporary_files merge_files;

    FastqReaderWrapper fastqReader;

    unsigned n_threads;

    bool verbose;

    std::mutex global_num_corrected_reads_mutex;

    std::size_t global_num_corrected_reads;

    std::size_t global_num_corrected_errors_step1_1;
    std::size_t global_num_corrected_errors_step1_2;
    std::size_t global_num_corrected_errors_step1_3;
    std::size_t global_num_corrected_errors_step2_1;
    std::size_t global_num_corrected_errors_step2_2;

    std::size_t global_num_corrected_substs;
    std::size_t global_num_corrected_ins;
    std::size_t global_num_corrected_dels;
    std::size_t global_num_corrected_pairs; // pairs of insertion and deletion in a single region
    std::size_t global_num_first_kmer_corrections;
    std::size_t global_num_first_kmer_successes;
    std::size_t global_num_single_deletions;
    std::size_t global_num_checked_with_long_kmer;
    std::size_t global_num_filtered_with_long_kmer;
    std::size_t global_num_instantly_accepted_with_long_kmer;

    std::size_t global_num_oriented_reads;

    std::size_t next_chunk;
    std::mutex next_chunk_mutex;

public:
    C_correct_errors(const C_arg& c_inst_args, const ReadFileData& _read_file_data, C_check_read& check_read, bool _accept_filtered_with_long_kmers) :
        c_log(std::cout),
        c_err(std::cerr),
        quality_score_offset(check_read.quality_score_offset),
        kmer_length(c_inst_args.kmer_length),
        long_kmer_length(c_inst_args.long_kmer_length),
        extend(c_inst_args.extend),
        accept_filtered_with_long_kmers(_accept_filtered_with_long_kmers),
        correct_indels(c_inst_args.correct_indels),
        b_change_headers_length_tag(c_inst_args.change_headers_length_tag),
        mark_corrected(c_inst_args.mark_corrected),
        read_file_data(_read_file_data),
        merge_files(_read_file_data),
        fastqReader(PART_SIZE, c_inst_args.n_threads * PART_BUFFERS_PER_THREAD),
        n_threads(c_inst_args.n_threads),
        verbose(c_inst_args.verbose),
        global_num_corrected_reads(0),
        global_num_corrected_errors_step1_1(0),
        global_num_corrected_errors_step1_2(0),
        global_num_corrected_errors_step1_3(0),
        global_num_corrected_errors_step2_1(0),
        global_num_corrected_errors_step2_2(0),
        global_num_corrected_substs(0),
        global_num_corrected_ins(0),
        global_num_corrected_dels(0),
        global_num_corrected_pairs(0),
        global_num_first_kmer_corrections(0),
        global_num_first_kmer_successes(0),
        global_num_single_deletions(0),
        global_num_checked_with_long_kmer(0),
        global_num_filtered_with_long_kmer(0),
        global_num_instantly_accepted_with_long_kmer(0),
        global_num_oriented_reads(0),
        next_chunk(0)
    {}

    C_correct_errors(const C_correct_errors&) = delete;

    void correct_errors_in_a_file(CKMCFile& kmc_file, CKMCFile* kmc_long_file);

    void operator()(CKMCFile& kmc_file, CKMCFile* kmc_long_file = nullptr);

private:
    template<bool CORRECT_INDELS>
    void change_headers_length_tag(std::string& header, std::size_t original_length, std::size_t new_length);
    void write_a_read_to_fastq(const std::string& header, const std::string& connector, const std::string& qualities, const std::string& sequence, FileReader& f_error_correction);

    template<bool CORRECT_INDELS>
    void correct_errors_in_a_chunk(CKMCFile& kmc_file, const std::size_t current_chunk_number, C_correct_read<CORRECT_INDELS>& correct_read, ReadsChunk& readsChunk);
};




//----------------------------------------------------------------------
// Reads lines from file and calls correction.
//----------------------------------------------------------------------

template<bool CORRECT_INDELS>
void C_correct_errors::correct_errors_in_a_chunk(CKMCFile& kmc_file, const std::size_t current_chunk_number, C_correct_read<CORRECT_INDELS>& correct_read, ReadsChunk& readsChunk) {
    std::string current_error_correction_info_file_name = read_file_data.getCorrectionInfoName(current_chunk_number);

    // open error correction information files
    FileReader f_error_correction;
    f_error_correction.setFileName(current_error_correction_info_file_name, read_file_data.type);

    // check error correction information files
    if (!f_error_correction.openFile(FileReader::WRITE)) {
        c_err << std::endl << "ERROR: Cannot create temporary file " << current_error_correction_info_file_name << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string original_sequence;
    std::string original_quality_score;

    // perform correction for every read in current chunk
    while (readsChunk.getLine(correct_read.header)) {

        // DNA sequence
        readsChunk.getLine(correct_read.sequence_modified);
        original_sequence = correct_read.sequence_modified;

        // "+"
        readsChunk.getLine(correct_read.connector);

        // quality score
        readsChunk.getLine(correct_read.quality_score);
        original_quality_score = correct_read.quality_score;

        bool too_many_errors(false);
        size_t current_read_length = correct_read.sequence_modified.length();
        bool read_modified = false;

        if (current_read_length > correct_read.kmer_length) {
            std::size_t number_of_Ns = correctionUtilities.replace_Ns(correct_read.sequence_modified);

            if (number_of_Ns >= current_read_length * MAX_N_RATIO) {
                too_many_errors = true;
            }

            //----------------------------------------------------------------------
            // correct errors in a read
            //----------------------------------------------------------------------

            // update num_corrected_reads
            if (!too_many_errors) {
                correct_read.prepare_corrector();

                std::size_t old_sum_reads = correct_read.num_corrected_errors_step1_1 + correct_read.num_corrected_errors_step1_2 + correct_read.num_corrected_errors_step1_3 + correct_read.num_corrected_errors_step2_1 + correct_read.num_corrected_errors_step2_2;

                correct_read.correct_errors_in_a_read();

                std::size_t new_sum_errors = correct_read.num_corrected_errors_step1_1 + correct_read.num_corrected_errors_step1_2 + correct_read.num_corrected_errors_step1_3 + correct_read.num_corrected_errors_step2_1 + correct_read.num_corrected_errors_step2_2;
                if ((new_sum_errors - old_sum_reads) > (current_read_length * MAX_ERROR_RATE)) {
                    too_many_errors = true;
                }
                else if ((new_sum_errors - old_sum_reads) > 0 || number_of_Ns > 0) {
                    correct_read.num_corrected_reads++;
                }
            }
        }

        if (b_change_headers_length_tag) {
            change_headers_length_tag<CORRECT_INDELS>(correct_read.header, original_sequence.length(), correct_read.sequence_modified.length());
            change_headers_length_tag<CORRECT_INDELS>(correct_read.connector, original_sequence.length(), correct_read.sequence_modified.length());
        }

        if (mark_corrected && read_modified) {
            correct_read.header += " corrected";
        }

        write_a_read_to_fastq(correct_read.header, correct_read.connector, too_many_errors ? original_quality_score : correct_read.quality_score, too_many_errors ? original_sequence : correct_read.sequence_modified, f_error_correction);
    }

    // close error correction information files
    f_error_correction.close();
}

#endif
