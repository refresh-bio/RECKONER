/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.2
 * 
 */

#include "correct_errors.hpp"
#include "FileReader.h"
#include "Log.h"
#include <kmc_api/kmc_file.h>
#include <algorithm>
#include <omp.h>
#include <cstring>

#define USE_THREADS
#define DONT_PRINT_COUNTER
#define COUNTER_FREQ 1

//----------------------------------------------------------------------
// Reads lines from file and paralelly calls correction.
//----------------------------------------------------------------------

void C_correct_errors::correct_errors_in_reads(const C_arg& c_inst_args, const std::string& error_correction_info_file_name, CKMCFile& kmc_file) {
    std::ofstream& f_log = Log::get_stream();
    std::size_t current_chunk = 0;
    std::size_t NN = 0;

#ifdef USE_THREADS
#pragma omp parallel num_threads(static_cast<int>(c_inst_args.n_threads))
#endif
    {
        std::size_t local_chunk = 0;
#ifdef USE_FILE_READER
        FileReader f_read;
#else
        std::ifstream f_read;
#endif

        std::size_t total_num_corrected_errors_step1_1(0);
        std::size_t total_num_corrected_errors_step1_2(0);
        std::size_t total_num_corrected_errors_step1_3(0);
        std::size_t total_num_corrected_errors_step2_1(0);
        std::size_t total_num_corrected_errors_step2_2(0);

        std::size_t num_corrected_errors_step1_1(0);
        std::size_t num_corrected_errors_step1_2(0);
        std::size_t num_corrected_errors_step1_3(0);
        std::size_t num_corrected_errors_step2_1(0);
        std::size_t num_corrected_errors_step2_2(0);

        std::size_t num_corrected_reads(0);

        C_correct_read correct_read(quality_score_offset, c_inst_args.kmer_length, c_inst_args.extend, max_read_length, read_file_name, num_corrected_errors_step1_1, num_corrected_errors_step1_2, num_corrected_errors_step1_3, num_corrected_errors_step2_1, num_corrected_errors_step2_2, kmc_file);

        const std::size_t n_chunks = chunks.size();

        // perform correction while some chunks are available
        while (true) {
            bool exit_loop = false;
#ifdef USE_THREADS
#pragma omp critical
#endif
            {
                if (current_chunk >= n_chunks) {
                    exit_loop = true;
                }
                local_chunk = current_chunk++;
            }

            if (exit_loop) {
                break;
            }

            // open error correction information files
            std::ofstream f_error_correction;

            f_error_correction.open(get_error_correction_info_file_name(error_correction_info_file_name, local_chunk).c_str(), std::ios::binary);

            // check error correction information files
            if (f_error_correction.is_open() == false) {
                std::cerr << std::endl << "ERROR: Cannot open " << get_error_correction_info_file_name(error_correction_info_file_name, local_chunk) << std::endl << std::endl;
                f_log << std::endl << "ERROR: Cannot open " << get_error_correction_info_file_name(error_correction_info_file_name, local_chunk) << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }

#ifdef USE_FILE_READER
            // open read files
            f_read.setFileName(read_file_name, read_file_type);

            // iterate reads
            if (f_read.openFile(FileReader::READ)) {
                f_read.seekPos(chunks[local_chunk]);

                // header
                f_read.getLine(correct_read.header);
#else
            // open read files
            f_read.open(read_file_name.c_str());

            // iterate reads
            if (f_read.is_open()) {
                f_read.seekg(chunks[local_chunk]);

                // header
                getline(f_read, correct_read.header);
#endif
                std::string& sequence_modification = correct_read.sequence_modification;


                std::size_t local_read = 0;

                // perform correction for every read in current chunk
                while ((local_chunk < (n_chunks - 1) && local_read < chunk_size) || (local_chunk == (n_chunks - 1) && local_read < last_chunk_size)) {
                    if (f_read.eof()) {
                        std::cerr << "ERROR: Erroneous chunk size (reads missing)." << std::endl;
                        f_log << "Erroneous chunk size (reads missing)." << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    num_corrected_errors_step1_1 = 0;
                    num_corrected_errors_step1_2 = 0;
                    num_corrected_errors_step1_3 = 0;
                    num_corrected_errors_step2_1 = 0;
                    num_corrected_errors_step2_2 = 0;

                    ++local_read;
#ifdef USE_FILE_READER
                    // DNA sequence
                    f_read.getLine(correct_read.sequence);
                    std::size_t current_read_length = correct_read.sequence.length();

                    // "+"
                    f_read.getLine(correct_read.connector);

                    // quality score
                    f_read.getLine(correct_read.quality_score);
#else
                    // DNA sequence
                    getline(f_read, correct_read.sequence);
                    std::size_t current_read_length = correct_read.sequence.length();

                    // "+"
                    getline(f_read, correct_read.connector);

                    // quality score
                    getline(f_read, correct_read.quality_score);
#endif
                    bool too_many_errors(false);
                    std::fill(sequence_modification.begin(), sequence_modification.begin() + current_read_length, '0');

                    if (current_read_length > correct_read.kmer_length) {
                        // change sequences to upper case
                        transform(correct_read.sequence.begin(), correct_read.sequence.end(), correct_read.sequence.begin(), ::toupper);

                        // count and substitute Ns other characters
                        std::size_t number_of_Ns = 0;
                        for (std::size_t it = 0; it < current_read_length; ++it) {
                            bool is_nucleotide = false;
                            for (std::size_t it_nucl = 0; it_nucl < NUM_NEOCLEOTIDE; ++it_nucl) {
                                if (correct_read.sequence[it] == NEOCLEOTIDE[it_nucl]) {
                                    is_nucleotide = true;
                                    break;
                                }
                            }
                            if (!is_nucleotide) {
                                ++number_of_Ns;
                                correct_read.sequence[it] = SUBST_CHAR;
                            }
                        }

                        ++NN;

                        //----------------------------------------------------------------------
                        // correct errors in a read
                        //----------------------------------------------------------------------
#ifdef USE_THREADS
#ifndef DONT_PRINT_COUNTER
#pragma omp critical
#endif
#endif
                        {
#ifndef DONT_PRINT_COUNTER
                            if (NN % COUNTER_FREQ == 0) {
                                std::cout << NN << "\r";
                                std::cout.flush();
                            }
#endif
                        }

                        // update num_corrected_reads
                        if (number_of_Ns < current_read_length * MAX_N_RATIO) {
                            correct_read.read_length = current_read_length;
                            correct_read.correct_errors_in_a_read_fastq();

                            std::size_t sum_errors = num_corrected_errors_step1_1 + num_corrected_errors_step1_2 + num_corrected_errors_step1_3 + num_corrected_errors_step2_1 + num_corrected_errors_step2_2;
                            if (sum_errors > (current_read_length * MAX_ERROR_RATE)) {
                                too_many_errors = true;
                            }
                            else if (sum_errors > 0) {
                                total_num_corrected_errors_step1_1 += num_corrected_errors_step1_1;
                                total_num_corrected_errors_step1_2 += num_corrected_errors_step1_2;
                                total_num_corrected_errors_step1_3 += num_corrected_errors_step1_3;
                                total_num_corrected_errors_step2_1 += num_corrected_errors_step2_1;
                                total_num_corrected_errors_step2_2 += num_corrected_errors_step2_2;

                                num_corrected_reads++;
                            }
                        }
                        else {
                            too_many_errors = true;
                        }
                    }

                    write_single_read(correct_read.sequence, sequence_modification, current_read_length, too_many_errors, f_error_correction);
#ifdef USE_FILE_READER
                    // header
                    f_read.getLine(correct_read.header);
#else
                    // header
                    getline(f_read, correct_read.header);
#endif
                }
            }
            else {
                std::cerr << std::endl << "ERROR: Cannot open " << read_file_name << std::endl << std::endl;
                f_log << std::endl << "ERROR: Cannot open " << read_file_name << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }

            // close read files
            f_read.close();

            // close error correction information files
            f_error_correction.close();
        }
#ifdef USE_THREADS
#pragma omp critical
#endif
        {
            global_num_corrected_reads += num_corrected_reads;
        }
    }
    std::cout << "     Number of corrected reads : " << std::setw(10) << global_num_corrected_reads << std::endl;
    std::cout << "     Correcting errors in reads: done" << std::endl << std::endl;

    f_log << "     Number of corrected reads : " << std::setw(10) << global_num_corrected_reads << std::endl;
    f_log << "     Correcting errors in reads: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// Saves corrected read into a temporary file.
//----------------------------------------------------------------------

void C_correct_errors::write_single_read(const std::string& correct_read_sequence, const std::string& sequence_modification, const std::size_t read_length, const bool too_many_errors, std::ofstream& f_error_correction) {
    // write remaining bps into the vectors
    const std::size_t residue_read(read_length % BPS_PER_BYTE);
    //----------------------------------------------------------------------
    // write error correction information to files
    //----------------------------------------------------------------------
    char buffer = ZERO;

    std::size_t it_mod;
    for (it_mod = 0; it_mod < read_length; it_mod++) {
        if (too_many_errors == true) {
            encode_correction_info(buffer, correct_read_sequence[it_mod], '0');
        }
        else {
            encode_correction_info(buffer, correct_read_sequence[it_mod], sequence_modification[it_mod]);
        }

        if ((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) {
            // update write_buffer*
            f_error_correction.write(&buffer, 1);

            // initialize buffer
            buffer &= ZERO;
        }
    }

    // read_length is a multiple of BPS_PER_BYTE
    // do nothing
    if (residue_read == 0) {
    } // read_length is not a multiple of BPS_PER_BYTE
    else {
        for (std::size_t it_fill = 0; it_fill < (BPS_PER_BYTE - residue_read); it_fill++) {
            buffer = buffer << 2;
        }

        // fill the last entry
        f_error_correction.write(&buffer, 1);
    }
}



//----------------------------------------------------------------------
// Saves read symbol in two bits.
//----------------------------------------------------------------------

inline void C_correct_errors::encode_correction_info(char& buffer, const char char_in_read, const char char_in_info) {
    buffer <<= 2;

    if (char_in_info != '0') {
        unsigned char read_id;
        unsigned char in_id;
        switch (char_in_info) {
            case 'A': in_id = 0;
                break;
            case 'C': in_id = 1;
                break;
            case 'G': in_id = 2;
                break;
            case 'T': in_id = 3;
                break;
            default: return;
        }

        switch (char_in_read) {
            case 'A': read_id = 0;
                break;
            case 'C': read_id = 1;
                break;
            case 'G': read_id = 2;
                break;
            case 'T': read_id = 3;
                break;
            default: return;
        }

        buffer |= (in_id - read_id + 4) % 4;
    }
}



//----------------------------------------------------------------------
// Extracts read symbol from two bit form.
//----------------------------------------------------------------------

inline char C_correct_errors::decode_correction_info(unsigned char first_bits, const char read) {
    std::ofstream& f_log = Log::get_stream();
    char symbols[] = {'A', 'C', 'G', 'T'};
    int symbol_id;

    if ((first_bits & (BIT8 | 0x40)) != first_bits) {
        std::cerr << std::endl << "ERROR: Illegal result in correction information file " << first_bits << std::endl << std::endl;
        f_log << std::endl << "ERROR: Illegal result in correction information file " << first_bits << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    switch (read) {
        case 'A': symbol_id = 0;
            break;
        case 'C': symbol_id = 1;
            break;
        case 'G': symbol_id = 2;
            break;
        case 'T': symbol_id = 3;
            break;
        default: return read;
    }

    first_bits >>= 6;
    return symbols[(symbol_id + first_bits) % 4];
}


//----------------------------------------------------------------------
// Merges correction data with input reads.
//----------------------------------------------------------------------

void C_correct_errors::write_corrected_reads(const C_arg& c_inst_args, const std::string& error_correction_info_file_name, const std::string& corrected_read_file_name, const FileReader::FileType corrected_read_file_type) {
    std::ofstream& f_log = Log::get_stream();

#ifdef USE_FILE_READER
    // open corrected read files
    FileReader f_corrected_read;
    FileReader f_read;
    f_corrected_read.setFileName(corrected_read_file_name, corrected_read_file_type);
    if (!f_corrected_read.openFile(FileReader::WRITE)) {
        std::cerr << std::endl << "ERROR: Cannot open " << corrected_read_file_name << std::endl << std::endl;
        f_log << std::endl << "ERROR: Cannot open " << corrected_read_file_name << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    // open read files
    f_read.setFileName(read_file_name, corrected_read_file_type);

    // process reads
    if (f_read.openFile(FileReader::READ)) {
#else
    // open corrected read files
    std::ofstream f_corrected_read;
    std::ifstream f_read;
    f_corrected_read.open(corrected_read_file_name.c_str());

    if (f_corrected_read.is_open() == false) {
        std::cerr << std::endl << "ERROR: Cannot open " << corrected_read_file_name << std::endl << std::endl;
        f_log << std::endl << "ERROR: Cannot open " << corrected_read_file_name << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    // open read files
    f_read.open(read_file_name.c_str());

    // process reads
    if (f_read.is_open()) {
#endif
        std::size_t next_chunk = 0;

        // correction information files
        std::ifstream f_error_correction_info;

        std::string line;
        std::string read;

#ifdef USE_FILE_READER
        // number of bytes for storing reads
        // std::size_t num_byte_per_read;
        // num_byte_per_read = (std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE));

        // std::string modification_buffer(num_byte_per_read, '0');
        while (f_read.getLine(line)) {
            // write the first head lines
            if (line.length() > 0) {
                line += "\n";
                write_line_successfully(f_corrected_read, line);
            }
            // DNA sequence
            f_read.getLine(read);
#else
        // header
        getline(f_read, line);

        // write the first head lines
        if (line.length() > 0) {
            f_corrected_read << line << "\n";
        }

        // number of bytes for storing reads
        // std::size_t num_byte_per_read;
        // num_byte_per_read = (std::size_t)(ceil((float)(read_length) / BPS_PER_BYTE));

        // std::string modification_buffer(num_byte_per_read, '0');

        while (!f_read.eof()) {
            // DNA sequence
            getline(f_read, read);
#endif
            const std::size_t current_read_length = read.length();

            // change sequences to upper case
            transform(read.begin(), read.end(), read.begin(), ::toupper);

            std::string original_read(read);

            // count and substitute Ns other characters
            std::size_t number_of_Ns = 0;
            for (std::size_t it = 0; it < current_read_length; ++it) {
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

            const bool too_many_Ns = number_of_Ns >= current_read_length * MAX_N_RATIO;

            // modified reads
            std::string read_modified(too_many_Ns ? original_read : read);

            // read modification information from files
            char buffer;

            if (!f_error_correction_info.get(buffer)) {
                f_error_correction_info.close();

                // there is no error correction information file
                if (next_chunk == chunks.size()) {
                    if (!f_read.eof()) {
                        std::cerr << std::endl << "ERROR: Unexpected end of error correction info file. " << std::endl;
                        f_log << std::endl << "ERROR: Unexpected end of error correction info file." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    // open error correction information files
                    f_error_correction_info.open(get_error_correction_info_file_name(error_correction_info_file_name, next_chunk).c_str(), std::ios::binary);

                    if (f_error_correction_info.is_open() == false) {
                        std::cerr << std::endl << "ERROR: Cannot open " << get_error_correction_info_file_name(error_correction_info_file_name, next_chunk) << std::endl << std::endl;
                        f_log << std::endl << "ERROR: Cannot open " << get_error_correction_info_file_name(error_correction_info_file_name, next_chunk) << std::endl << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    if (!f_error_correction_info.get(buffer)) {
                        std::cerr << std::endl << "ERROR: The file " << get_error_correction_info_file_name(error_correction_info_file_name, next_chunk) << " is empty" << std::endl << std::endl;
                        f_log << std::endl << "ERROR: The file " << get_error_correction_info_file_name(error_correction_info_file_name, next_chunk) << " is empty" << std::endl << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }

                // remove previous file
                if (next_chunk != 0) {
                    remove_error_correction_info_file(error_correction_info_file_name, next_chunk - 1);
                }
                ++next_chunk;

            }

            // make new reads
            std::size_t it_mod;
            for (it_mod = 0; it_mod < current_read_length; it_mod++) {
                unsigned char first_bits = buffer & (BIT8 | 0x40);
                buffer <<= 2;

                if (first_bits != 0 && !too_many_Ns) {
                    read_modified[it_mod] = decode_correction_info(first_bits, read[it_mod]);
                }

                // increment indexes
                if (((it_mod % BPS_PER_BYTE) == (BPS_PER_BYTE - 1)) && (it_mod != (current_read_length - 1))) {
                    if (!f_error_correction_info.get(buffer)) {
                        std::cerr << std::endl << "ERROR: The size of " << error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                        f_log << std::endl << "ERROR: The size of " << error_correction_info_file_name << " is wrong" << std::endl << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
#ifdef USE_FILE_READER
            // write new reads to output files
            read_modified += "\n";
            write_line_successfully(f_corrected_read, read_modified);

            // "+"
            f_read.getLine(line);
            line += "\n";
            write_line_successfully(f_corrected_read, line);

            // quality score
            f_read.getLine(line);
            line += "\n";
            write_line_successfully(f_corrected_read, line);
#else
            // write new reads to output files
            f_corrected_read << read_modified << "\n";

            // "+"
            getline(f_read, line);
            f_corrected_read << line << "\n";

            // quality score
            getline(f_read, line);
            f_corrected_read << line << "\n";

            // header
            getline(f_read, line);

            if (line.length() > 0) {
                f_corrected_read << line << "\n";
            }
#endif
        }

        f_error_correction_info.close();

        // remove last file
        remove_error_correction_info_file(error_correction_info_file_name, next_chunk - 1);
    }

    // close read files
    f_read.close();

    // close corrected reads
    f_corrected_read.close();

    std::cout << "     Writing corrected reads into files: done" << std::endl << std::endl;
    f_log << "     Writing corrected reads into files: done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// Removes unnecessary temporary file.
//----------------------------------------------------------------------

void C_correct_errors::remove_error_correction_info_file(const std::string& error_correction_info_file_name, const std::size_t file_index) {
    remove(get_error_correction_info_file_name(error_correction_info_file_name, file_index).c_str());
}




//----------------------------------------------------------------------
// Writes line to file and checks the operation success.
//----------------------------------------------------------------------

inline void C_correct_errors::write_line_successfully(FileReader& output_file, const std::string& line) {
    if (!output_file.putString(line)) {
        std::ofstream& f_log = Log::get_stream();
        std::string output_file_name;
        output_file.getFileName(output_file_name);

        std::cerr << std::endl << "ERROR: Cannot write to file " << output_file_name << ". Probably disk is full." << std::endl << std::endl;
        f_log << std::endl << "ERROR: Cannot write to file " << output_file_name << ". Probably disk is full." << std::endl << std::endl;

        exit(EXIT_FAILURE);
    }
}




//----------------------------------------------------------------------
// Generates temporary file name.
//----------------------------------------------------------------------

std::string C_correct_errors::get_error_correction_info_file_name(const std::string& error_correction_info_file_name, const std::size_t file_number) const {
    std::ostringstream error_correction_info_file_name_stream;
    error_correction_info_file_name_stream << error_correction_info_file_name;
    error_correction_info_file_name_stream << file_number;

    return error_correction_info_file_name_stream.str();
}
