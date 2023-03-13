/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#ifndef _PARSE_ARGS_H
#define _PARSE_ARGS_H

#include "define.hpp"
#include "FileReader.h"
#include <vector>



struct ReadFileData {
    std::string input_name;
    FileReader::FileType type;
    std::string extension;
    std::string raw_name;
    std::string correction_info_name;
    std::string output_name;

    ReadFileData(const std::string& _input_name);

    void addPrefix(const std::string& prefix);
    std::string getCorrectionInfoName(const std::size_t fileNumber) const;
};



//----------------------------------------------------------------------
// C_arg
//----------------------------------------------------------------------

class C_arg {
public:
    // variables
    std::vector<ReadFileData> read_files_data;

    std::string kmc_determine_params_database_name;
    std::string kmc_database_name;
    std::string kmc_long_database_name;
    std::string kmc_filtered_database_name;
    std::string kmc_list_file_name;
    std::string log_file_name;

    std::string prefix;

    bool is_kmer_length_user_defined;
    mutable std::size_t kmer_length;
    mutable std::size_t long_kmer_length;

    bool is_genome_size_user_defined;
    mutable std::size_t genome_size;

    std::size_t extend;

    unsigned n_threads;
    std::size_t kmc_memory;
    bool kmc_ram;

    bool reuse_kmc_db;
    bool use_long_kmer;
    double long_kmer_ratio;
    bool accept_filtered_with_long_kmers;

    bool correct_indels;
    bool change_headers_length_tag;

    bool mark_corrected;
    bool verbose;

    // constructors

    explicit C_arg(int argc, char** argv);

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
