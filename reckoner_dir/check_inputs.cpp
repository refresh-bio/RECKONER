/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.1.1
 * 
 */

#include "check_inputs.hpp"
#include "FileReader.h"
#include "define.hpp"
#include "Log.h"
#include <string>
#include <iostream>



//----------------------------------------------------------------------
// Reads and analyzes entire input file.
//----------------------------------------------------------------------

void C_check_read::check_read_file(const C_arg& c_inst_args, const std::string& read_file_name, const FileReader::FileType read_file_type) {
    std::ofstream& f_log = Log::get_stream();

    FileReader f_read;
    f_read.setFileName(read_file_name, read_file_type);

    if (!f_read.openFile(FileReader::READ)) {
        std::cerr << "ERROR: Cannot open " << read_file_name << " for checking reads." << std::endl;
        f_log << "ERROR: Cannot open " << read_file_name << " for checking reads." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    bool quality_score_offset_found = false;

    // header
    while (!quality_score_offset_found && f_read.getLine(line)) {
        // sequence
        f_read.getLine(line);
        // connector
        f_read.getLine(line);
        // quality score
        f_read.getLine(line);

        for (std::size_t it = 0; it < line.length(); it++) {
            if ((int)line[it] > MAX_64_SCORE - MAX_DETECT_SCORE_DIFF) {
                quality_score_offset = PHRED64;
                quality_score_offset_found = true;
                break;
            }

            if ((int)line[it] < MIN_33_SCORE + MAX_DETECT_SCORE_DIFF) {
                quality_score_offset = PHRED33;
                quality_score_offset_found = true;
                break;
            }
        }
    }
    std::cout << "     Quality score offset          : " << quality_score_offset << std::endl;
    f_log << "     Quality score offset          : " << quality_score_offset << std::endl;
}
