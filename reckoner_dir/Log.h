/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.1
 * 
 */

#ifndef LOG_H
#define LOG_H

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



class C_log {
private:
    std::ostream& stream1;
    std::ostream& stream2;

public:
    C_log(std::ostream& _stream1) : stream1(_stream1), stream2(Log::get_stream()) {}

    template<typename Type>
    C_log& operator<< (const Type& data) {
        stream1 << data;
        stream2 << data;

        return *this;
    }

    //for std::endl support
    typedef std::basic_ostream<char, std::char_traits<char> > StdStream;
    C_log& operator<<(StdStream& (*fp)(StdStream&)) {
        fp(stream1);
        fp(stream2);

        return *this;
    }
};




#endif /* LOG_H */

