/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.1.1
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
// Generates temporary file name.
//----------------------------------------------------------------------

std::string C_correction_utilities::get_error_correction_info_file_name(const std::string& error_correction_info_file_name, const std::string& output_file_name, const std::size_t file_number) const {
    if (file_number == 0) {
        return output_file_name;
    }

    std::ostringstream error_correction_info_file_name_stream;
    error_correction_info_file_name_stream << error_correction_info_file_name;
    error_correction_info_file_name_stream << file_number;

    return error_correction_info_file_name_stream.str();
}



//----------------------------------------------------------------------
// Replaces Ns in a read
//----------------------------------------------------------------------

bool C_correction_utilities::replace_Ns(std::string& read) const {
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

    const bool too_many_Ns = number_of_Ns >= read.length() * MAX_N_RATIO;
    return too_many_Ns;
}



//----------------------------------------------------------------------
// If number of Ns is below the threshold, removes Ns
//----------------------------------------------------------------------

bool C_correction_utilities::replace_Ns_revert(std::string& read) const {
    std::string originalRead(read);

    const bool too_many_Ns = replace_Ns(read);

    if (too_many_Ns) {
        read.swap(originalRead);
    }

    return too_many_Ns;
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
            output_file.open(output_file_name, std::ios_base::binary | std::ios_base::app);

            if (!output_file.is_open()) {
                std::cerr << std::endl << "ERROR: Cannot create output file." << std::endl << std::endl;
                f_log << std::endl << "ERROR: Cannot create output file." << std::endl << std::endl;
            }
        }
        // Append one file.
        else {
            std::string current_correction_info_file_name = correction_utilities.get_error_correction_info_file_name(correction_info_file_name, output_file_name, next_chunk);

            std::ifstream input_file(current_correction_info_file_name, std::ios_base::binary);
            if (!input_file.is_open()) {
                std::cerr << std::endl << "ERROR: Cannot open temporary file " << current_correction_info_file_name << std::endl << std::endl;
                f_log << std::endl << "ERROR: Cannot open temporary file " << current_correction_info_file_name << std::endl << std::endl;
            }

            output_file << input_file.rdbuf();

            if (!output_file) {
                std::cerr << std::endl << "ERROR: Cannot write to file " << output_file_name << ". Probably disk is full." << std::endl << std::endl;
                f_log << std::endl << "ERROR: Cannot write to file " << output_file_name << ". Probably disk is full." << std::endl << std::endl;

                exit(EXIT_FAILURE);
            }

            input_file.close();
            remove_error_correction_info_file(correction_info_file_name, output_file_name, next_chunk);
        }
        ++next_chunk;
    }
}



//----------------------------------------------------------------------
// Creates correction threads.
//----------------------------------------------------------------------

void C_correct_errors::correct_errors_in_reads(CKMCFile& kmc_file) {
    if (!fastqReader.openFile(input_file_name, file_type)) {
        std::cerr << std::endl << "ERROR: Cannot open input file" << std::endl << std::endl;
        f_log << std::endl << "ERROR: Cannot open input file" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    std::thread reader_thread(std::ref(fastqReader));
    std::thread merger_thread(std::ref(merge_files));

    std::vector<std::thread> threads(n_threads);

    for (size_t it = 0; it < threads.size(); ++it) {
        threads[it] = std::thread(std::ref(*this), std::ref(kmc_file));
    }

    for (size_t it = 0; it < threads.size(); ++it) {
        threads[it].join();
    }

    merge_files.signal_finish();
    merger_thread.join();
    reader_thread.join();

    std::cout << "     Number of corrected reads     : " << global_num_corrected_reads << std::endl << std::endl;
    f_log << "     Number of corrected reads     : " << global_num_corrected_reads << std::endl << std::endl;
}



//----------------------------------------------------------------------
// Correction thread body.
//----------------------------------------------------------------------

void C_correct_errors::operator()(CKMCFile& kmc_file) {
    C_correct_read correct_read(quality_score_offset, kmer_length, extend, input_file_name, kmc_file);

    // perform correction while some chunks are available
    while (true) {
        ReadsChunk readsChunk;
        if (!fastqReader.getChunk(readsChunk)) {
            break;
        }

        correct_errors_in_reads_single_chunk(error_correction_info_file_name, kmc_file, readsChunk.getChunkNo(), correct_read, readsChunk);

        merge_files.signal_file_availability(readsChunk.getChunkNo());
    }

    { // for lock_guard
        std::lock_guard<std::mutex> global_num_corrected_reads_guard(global_num_corrected_reads_mutex);
        global_num_corrected_reads += correct_read.num_corrected_reads;
    }
}




//----------------------------------------------------------------------
// Saves corrected read into a temporary file.
//----------------------------------------------------------------------

void C_correct_errors::write_single_read_to_fastq(const std::string& header, const std::string& connector, const std::string& qualities, std::string& original_sequence, const std::string& correct_read_sequence, const std::string& sequence_modification, const std::size_t read_length, const bool too_many_errors, FileReader& f_error_correction) {
    f_error_correction.putString(header);
    f_error_correction.putString("\n");

    if (correctionUtilities.replace_Ns_revert(original_sequence)) {
        f_error_correction.putString(original_sequence);
        f_error_correction.putString("\n");
    }
    else if (too_many_errors) {
        f_error_correction.putString(correct_read_sequence);
        f_error_correction.putString("\n");
    }
    else {
        std::string corrected_sequence = correct_read_sequence;
        for (unsigned i = 0; i < correct_read_sequence.length(); ++i) {
            if (sequence_modification[i] != '0') {
                corrected_sequence[i] = sequence_modification[i];
            }
        }
        f_error_correction.putString(corrected_sequence);
        f_error_correction.putString("\n");
    }

    f_error_correction.putString(connector);
    f_error_correction.putString("\n");
    f_error_correction.putString(qualities);
    f_error_correction.putString("\n");
}



//----------------------------------------------------------------------
// Reads lines from file and calls correction.
//----------------------------------------------------------------------

void C_correct_errors::correct_errors_in_reads_single_chunk(const std::string& error_correction_info_file_name, CKMCFile& kmc_file, const std::size_t current_chunk_number, C_correct_read& correct_read, ReadsChunk& readsChunk) {
    std::string current_error_correction_info_file_name = correctionUtilities.get_error_correction_info_file_name(error_correction_info_file_name, output_file_name, current_chunk_number);

    // open error correction information files
    FileReader f_error_correction;
    f_error_correction.setFileName(current_error_correction_info_file_name, file_type);

    // check error correction information files
    if (!f_error_correction.openFile(FileReader::WRITE)) {
        std::cerr << std::endl << "ERROR: Cannot create temporary file " << current_error_correction_info_file_name << std::endl << std::endl;
        f_log << std::endl << "ERROR: Cannot create temporary file " << current_error_correction_info_file_name << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string original_sequence;

    // perform correction for every read in current chunk
    while (readsChunk.getLine(correct_read.header)) {

        // DNA sequence
        readsChunk.getLine(correct_read.sequence);
        original_sequence = correct_read.sequence;

        size_t current_read_length = correct_read.sequence.length();
        correct_read.set_read_length(current_read_length);

        // "+"
        readsChunk.getLine(correct_read.connector);

        // quality score
        readsChunk.getLine(correct_read.quality_score);

        bool too_many_errors(false);

        if (current_read_length > correct_read.kmer_length) {
            too_many_errors = correctionUtilities.replace_Ns(correct_read.sequence);

            //----------------------------------------------------------------------
            // correct errors in a read
            //----------------------------------------------------------------------

            // update num_corrected_reads
            if (!too_many_errors) {
                correct_read.prepare_corrector();
                correct_read.correct_errors_in_a_read_fastq();

                std::size_t sum_errors = correct_read.num_corrected_errors_step1_1 + correct_read.num_corrected_errors_step1_2 + correct_read.num_corrected_errors_step1_3 + correct_read.num_corrected_errors_step2_1 + correct_read.num_corrected_errors_step2_2;
                if (sum_errors > (current_read_length * MAX_ERROR_RATE)) {
                    too_many_errors = true;
                }
                else if (sum_errors > 0) {
                    correct_read.num_corrected_reads++;
                }
            }
        }

        write_single_read_to_fastq(correct_read.header, correct_read.connector, correct_read.quality_score, original_sequence, correct_read.sequence, correct_read.sequence_modification, current_read_length, too_many_errors, f_error_correction);
    }

    // close error correction information files
    f_error_correction.close();
}
