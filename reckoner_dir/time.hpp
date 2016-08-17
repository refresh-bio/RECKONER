/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.2
 * 
 */

#ifndef _TIME_LOCAL_H
#define _TIME_LOCAL_H



#include "define.hpp"
#include <vector>



//----------------------------------------------------------------------
// C_time
//----------------------------------------------------------------------

class C_time {
public:
    // variables
    // parse arguments
    std::string start_parse_args;
    std::string end_parse_args;

    // check read files
    std::vector<std::string> vector_start_check_read_file;
    std::vector<std::string> vector_end_check_read_file;

    // correct errors in reads
    std::vector<std::string> vector_start_correct_errors_in_reads;
    std::vector<std::string> vector_end_correct_errors_in_reads;

    // write corrected reads
    std::vector<std::string> vector_start_write_corrected_reads;
    std::vector<std::string> vector_end_write_corrected_reads;

    // constructors

    C_time() {
    };
};



#endif
