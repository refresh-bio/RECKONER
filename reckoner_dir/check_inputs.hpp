/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.2
 * 
 */

#ifndef _PARSE_READ_H
#define _PARSE_READ_H

#include "parse_args.hpp"
#include "time.hpp"



//----------------------------------------------------------------------
// C_check_read
//----------------------------------------------------------------------

class C_check_read {
public:
    // variables
    std::size_t quality_score_offset;

    // constructor

    C_check_read() :
        quality_score_offset(PHRED33) {
    };

    // functions
    void check_read_file(const C_arg& c_inst_args, const std::string& read_file_name, const FileReader::FileType read_file_type);
};



#endif
