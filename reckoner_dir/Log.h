/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#ifndef LOG_H
#define	LOG_H

#include <fstream>
#include <string>

//----------------------------------------------------------------------
// Class methods are not thread safe.
//----------------------------------------------------------------------

class Log {
    Log() = delete;
    static std::ofstream f_log;
public:
    static void open_log_file(const std::string& file_name);
    static std::ofstream& get_stream();
};

#endif	/* LOG_H */

