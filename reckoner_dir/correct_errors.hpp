/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.2
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
    std::string get_error_correction_info_file_name(const std::string& error_correction_info_file_name, const std::string& output_file_name, const std::size_t file_number) const;

    bool replace_Ns(std::string& read) const;
    bool replace_Ns_revert(std::string& read) const;
};



class C_merge_temporary_files {
private:
    std::ofstream& f_log;

    C_correction_utilities correction_utilities;

    const std::string output_file_name;
    const std::string correction_info_file_name;

    std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> waiting_chunks;
    std::mutex waiting_chunks_mutex;
    std::condition_variable available_chunks;

    bool merging_finished;

    void remove_error_correction_info_file(const std::string& error_correction_info_file_name, const std::string& corrected_read_file_name, const std::size_t file_index) const {
        remove(correction_utilities.get_error_correction_info_file_name(error_correction_info_file_name, corrected_read_file_name, file_index).c_str());
    }

public:
    C_merge_temporary_files(const std::string& _output_file_name, const std::string& _correction_info_file_name) :
        f_log(Log::get_stream()),
        output_file_name(_output_file_name),
        correction_info_file_name(_correction_info_file_name),
        merging_finished(false)
    {}

    void signal_file_availability(size_t file_number);
    void signal_finish();

    void operator()();
};




class C_correct_errors {
private:
    std::ofstream& f_log;

    C_correction_utilities correctionUtilities;

    std::size_t quality_score_offset;
    std::size_t kmer_length;
    std::size_t extend;

    std::string error_correction_info_file_name;

    std::string input_file_name;
    std::string output_file_name;
    FileReader::FileType file_type;

    C_merge_temporary_files merge_files;

    FastqReaderWrapper fastqReader;

    unsigned n_threads;

    std::size_t global_num_corrected_reads;
    std::mutex global_num_corrected_reads_mutex;

    std::size_t next_chunk;
    std::mutex next_chunk_mutex;

public:
    C_correct_errors(const C_arg& c_inst_args, const std::string& _error_correction_info_file_name, const std::string& _input_file_name, const std::string& _output_file_name, FileReader::FileType _file_type, C_check_read& check_read) :
        f_log(Log::get_stream()),
        quality_score_offset(check_read.quality_score_offset),
        kmer_length(c_inst_args.kmer_length),
        extend(c_inst_args.extend),
        error_correction_info_file_name(_error_correction_info_file_name),
        input_file_name(_input_file_name),
        output_file_name(_output_file_name),
        file_type(_file_type),
        merge_files(_output_file_name, _error_correction_info_file_name),
        fastqReader(PART_SIZE, c_inst_args.n_threads * PART_BUFFERS_PER_THREAD),
        n_threads(c_inst_args.n_threads),
        global_num_corrected_reads(0),
        next_chunk(0)
    {}

    C_correct_errors(const C_correct_errors&) = delete;

    void correct_errors_in_reads(CKMCFile& kmc_file);

    void operator()(CKMCFile& kmc_file);

private:
    void write_single_read_to_fastq(const std::string& header, const std::string& connector, const std::string& qualities, std::string& original_sequence, const std::string& correct_read_sequence, const std::string& sequence_modification, const std::size_t read_length, const bool too_many_errors, FileReader& f_error_correction);

    void correct_errors_in_reads_single_chunk(const std::string& error_correction_info_file_name, CKMCFile& kmc_file, const std::size_t current_chunk_number, C_correct_read& correct_read, ReadsChunk& readsChunk);
};



#endif
