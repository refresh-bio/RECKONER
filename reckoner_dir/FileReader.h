/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#ifndef _FILE_READER_H
#define _FILE_READER_H

#include "Log.h"
#include <zlib.h>
#include <cstdlib>
#include <string>
#include <iostream>



#define BUFFER_SIZE (8 << 20)
#define GZIP_BUFFER_SIZE (8 << 20)

class FileReader {
    C_log c_err;

public:

    enum FileType {
        NOT_DEFINED, RAW, GZIP
    };

    enum OpenMode {
        READ, WRITE
    };

    bool getLine(char* outBuff);
    bool getLine(std::string& outString);
    bool putString(const std::string& buff);

    void setFileName(const std::string& _fileName, FileType _fileType);
    bool openFile(OpenMode openMode);
    void close();

    void getFileName(std::string& _fileName) {
        _fileName = fileName;
    }

    long tellPos();
    long seekPos(long pos);
    bool eof() const;

    FileReader() : c_err(std::cerr), file(NULL), gz_file(NULL), fileType(NOT_DEFINED), leaveSymbols(0), currentPos(0), fileBytesRead(0), eofReached(false), windowsNewLine(false) {
        buffer = new char[BUFFER_SIZE + 1];
    }

    ~FileReader() {
        close();
        delete[] buffer;
    }

private:
    FILE* file;
    gzFile_s* gz_file;

    std::string fileName;
    FileType fileType;

    std::size_t leaveSymbols;
    std::size_t currentPos;
    std::size_t fileBytesRead;

    bool eofReached;
    bool windowsNewLine;

    char* buffer;

    bool getEolPos(char* buff, std::size_t maxCharacters, std::size_t& outEolPos);
    std::size_t fileRead(char* outBuff, std::size_t elemSize, std::size_t elemCount, FILE* file);
};

#endif
