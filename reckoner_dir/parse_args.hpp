/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.2
 * 
 */

#ifndef _PARSE_ARGS_H
#define _PARSE_ARGS_H

#include "define.hpp"
#include "FileReader.h"
#include <vector>



//----------------------------------------------------------------------
// C_arg
//----------------------------------------------------------------------

class C_arg {
public:
    // variables
    std::vector<std::string> read_files_names;
    std::vector<FileReader::FileType> read_files_types;

    std::string kmc_determine_params_database_name;
    std::string kmc_database_name;
    std::string kmc_filtered_database_name;
    std::string kmc_list_file_name;
    std::string log_file_name;

    std::vector<std::string> error_correction_info_files_names;
    std::vector<std::string> corrected_read_files_names;

    std::string prefix;

    bool is_kmer_length_user_defined;
    mutable std::size_t kmer_length;

    bool is_genome_size_user_defined;
    mutable std::size_t genome_size;

    std::size_t extend;

    unsigned n_threads;
    std::size_t kmc_memory;

    // constructors

    explicit C_arg(int argc, char** argv);

    void extract_name_and_extension(std::string input_file_name, std::string& output_file_name, std::string& output_full_extension, FileReader::FileType& output_file_type);

private:
    // variables
    int num_args;

    char** args;

    // functions
    void read_args();
    void print_usage();
    void params_at_the_beginning();
};



#endif
