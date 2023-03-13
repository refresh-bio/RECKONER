/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#ifndef _PARSE_READ_H
#define _PARSE_READ_H

#include "parse_args.hpp"
#include "time.hpp"
#include "Log.h"



//----------------------------------------------------------------------
// C_check_read
//----------------------------------------------------------------------

class C_check_read {
    C_log c_log;
    C_log c_err;

public:
    // variables
    std::size_t quality_score_offset;

    // constructor

    C_check_read() :
        c_log(std::cout),
        c_err(std::cerr),
        quality_score_offset(PHRED33) {
    };

    // functions
    void check_read_file(const ReadFileData& read_file_data);
};



#endif
