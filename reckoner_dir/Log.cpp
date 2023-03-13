/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#include "Log.h"
#include <iostream>
#include <cstdlib>



//----------------------------------------------------------------------
// Opens file to save correction logs.
//----------------------------------------------------------------------

void Log::open_log_file(const std::string& file_name) {
    if (f_log.is_open()) {
        std::cerr << "ERROR: Cannot open log file - file has been already opened." << std::endl;
        exit(EXIT_FAILURE);
    }

    f_log.open(file_name.c_str());

    if (!f_log.is_open()) {
        std::cerr << "ERROR: Cannot create log file. Probably RECKONER has no write permission in the current folder" << std::endl;
        std::cerr << " or disk is completely full." << std::endl;
        exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------
// Returns the log file output stream.
//----------------------------------------------------------------------

std::ofstream& Log::get_stream() {
    if (!f_log.is_open()) {
        std::cerr << "ERROR: Log file has not been opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    return f_log;
}

std::ofstream Log::f_log;
