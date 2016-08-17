/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.2
 * 
 */

#ifndef _CORRECT_H
#define _CORRECT_H



#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "correct_read.hpp"
#include "FileReader.h"
#include <string>



//----------------------------------------------------------------------
// C_correct_errors
//----------------------------------------------------------------------

class C_correct_errors {
public:
    // variables
    std::size_t quality_score_offset;
    std::size_t num_reads;
    /* const */ std::size_t max_read_length; // leak of const - C_correct_erros have to be CopyAssignable

    std::vector<std::streampos> chunks;
    std::size_t chunk_size;
    std::size_t last_chunk_size;

    // constructors

    C_correct_errors(const C_arg& c_inst_args, C_check_read& check_read, const std::string& _read_file_name, const FileReader::FileType _read_file_type) :
        quality_score_offset(check_read.quality_score_offset),
        num_reads(check_read.num_reads),
        max_read_length(check_read.max_read_length),
        chunks(check_read.chunks),
        chunk_size(check_read.chunk_size),
        last_chunk_size(check_read.last_chunk_size),
        read_file_name(_read_file_name),
        read_file_type(_read_file_type),
        global_num_corrected_reads(0) {
    };

    // functions
    void correct_errors_in_reads(const C_arg& c_inst_args, const std::string& error_correction_info_file_name, CKMCFile& kmc_file);

    void write_corrected_reads(const C_arg& c_inst_args, const std::string& error_correction_info_file_name, const std::string& corrected_read_file_name, const FileReader::FileType corrected_read_file_type);

private:
    // variables
    std::string read_file_name;
    FileReader::FileType read_file_type;

    std::size_t global_num_corrected_reads;

    //functions
    void encode_correction_info(char& buffer, const char char_in_read, const char char_in_info);

    void write_single_read(const std::string& correct_read_sequence, const std::string& sequence_modification, const std::size_t read_length, const bool too_many_errors, std::ofstream& f_error_correction);

    char decode_correction_info(unsigned char first_bits, const char read);

    void remove_error_correction_info_file(const std::string& error_correction_info_file_name, const std::size_t file_index);

    void write_line_successfully(FileReader& output_file, const std::string& line);

    std::string decode_a_byte(char& in, const std::size_t& num_empty_characters);

    std::string get_error_correction_info_file_name(const std::string& error_correction_info_file_name, const std::size_t file_number) const;
};



#endif
