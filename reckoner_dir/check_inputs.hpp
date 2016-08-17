/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.2
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
    std::size_t num_reads;
    std::size_t read_length;
    std::size_t max_read_length;
    std::size_t init_num_elements;
    std::size_t interval;
    std::size_t min_num_processed_reads;
    std::size_t quality_score_offset;

    int max_quality_score;
    int min_quality_score;

    std::vector<std::streampos> chunks;
    std::size_t chunk_size;
    std::size_t last_chunk_size;

    // constructor

    C_check_read() :
        num_reads(0),
        max_read_length(0),
        init_num_elements(0),
        interval(0),
        min_num_processed_reads(0),
        quality_score_offset(0),
        max_quality_score(-1000),
        min_quality_score(1000),
        chunk_size (DEFAULT_CHUNK_SIZE),
        last_chunk_size(0) {
    };

    // functions
    void check_read_file(const C_arg& c_inst_args, const std::string& read_file_name, const FileReader::FileType read_file_type);
};



#endif
