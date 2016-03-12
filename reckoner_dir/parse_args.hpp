/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#ifndef _PARSE_ARGS_H
#define _PARSE_ARGS_H



#include "define.hpp"
#include "time.hpp"
#include "FileReader.h"
#include <vector>
#include <omp.h>



//----------------------------------------------------------------------
// C_arg
//----------------------------------------------------------------------

class C_arg {
public:
    // variables
    std::vector<std::string> read_files_names;
    std::vector<FileReader::FileType> read_files_types;
    std::string kmc_database_name;
    std::string log_file_name;
    std::string prefix;
    std::vector<std::string> error_correction_info_files_names; // nazwa pliku tymczasowego
    std::vector<std::string> corrected_read_files_names; // nazwa pliku wyjściowego

    std::size_t kmer_length;
    std::size_t extend;

    std::size_t n_threads;

    bool paired_read;
    bool nowrite;

    // constructors

    explicit C_arg(int argc, char** argv, C_time& c_inst_time) :
    kmer_length(0),
    extend(MAX_EXTENSION),
    n_threads(omp_get_max_threads()),
    paired_read(true),
    nowrite(false),
    num_args(argc),
    args(argv) {
        time_t rawtime;
        time(&rawtime);
        c_inst_time.start_parse_args = asctime(localtime(&rawtime));

        read_args();

        time(&rawtime);
        c_inst_time.end_parse_args = asctime(localtime(&rawtime));
    };

private:
    // variables
    int num_args;

    char** args;

    std::ofstream f_log;

    // functions
    void read_args();
    void print_usage();
};



#endif
