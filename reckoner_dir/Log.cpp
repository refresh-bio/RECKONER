/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
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
        std::cout << "ERROR: Cannot open log file - file has been already opened." << std::endl;
        exit(EXIT_FAILURE);
    }

    f_log.open(file_name.c_str());

    if (!f_log.is_open()) {
        std::cout << "ERROR: Cannot open log file." << std::endl;
        exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------
// Returns the log file output stream.
//----------------------------------------------------------------------

std::ofstream& Log::get_stream() {
    if (!f_log.is_open()) {
        std::cout << "ERROR: Log file has not been opened." << std::endl;
        exit(EXIT_FAILURE);
    }
    return f_log;
}

std::ofstream Log::f_log;
