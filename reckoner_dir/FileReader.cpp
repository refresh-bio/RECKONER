/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.1.1
 * 
 */

#define _CRT_SECURE_NO_DEPRECATE

#include "FileReader.h"
#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cassert>



//----------------------------------------------------------------------
// Finds the EOL in a buffer.
//----------------------------------------------------------------------

bool FileReader::getEolPos(char* buff, std::size_t maxCharacters, std::size_t& outEolPos) {
    leaveSymbols = 0;

    for (std::size_t i = 0; i < maxCharacters; ++i) {
        if (buff[i] == '\r') {
            outEolPos = i;

            if (i == maxCharacters - 1) {
                windowsNewLine = true;
            }
            else if (buff[i + 1] == '\n') {
                leaveSymbols = 1;
            }
            return true;
        }
        else if (buff[i] == '\n') {
            outEolPos = i;
            return true;
        }
    }

    outEolPos = maxCharacters;
    return false;
}

//----------------------------------------------------------------------
// Reads a line from opened file.
// outBuff must be of size >= BUFFER_SIZE
//----------------------------------------------------------------------

bool FileReader::getLine(char* outBuff) {
    if (eofReached) {
        return false;
    }

    currentPos += leaveSymbols;

    buffer[BUFFER_SIZE] = 0;

    std::size_t out;
    if (getEolPos(buffer + currentPos, fileBytesRead - currentPos, out)) {
        if (out + currentPos >= BUFFER_SIZE) {
            std::cerr << "Wrong file format." << std::endl;
            exit(EXIT_FAILURE);
        }
        memcpy(outBuff, buffer + currentPos, out);
        outBuff[out] = 0;
        currentPos += out + 1;
        return true;
    }
    else {
        // reached end of buffer
        std::size_t copyBytes = fileBytesRead - currentPos;
        memcpy(outBuff, buffer + currentPos, copyBytes);
        fileBytesRead = fileRead(buffer, sizeof (char), BUFFER_SIZE, file);

        currentPos = 0;

        if (copyBytes == 0 && windowsNewLine) {
            if (buffer[0] == '\n') {
                ++currentPos; // like leaveSymbols == 1
            }
            else {
                assert(false);
            }
            windowsNewLine = false;
        }

        const bool found = getEolPos(buffer + currentPos, fileBytesRead - currentPos, out);

        if (out + copyBytes >= BUFFER_SIZE - 1) {
            std::cerr << "Wrong file format." << std::endl;
            exit(EXIT_FAILURE);
        }

        if (!found) {
            eofReached = true;

            if (out + copyBytes == 0) {
                return false;
            }
        }

        memcpy(outBuff + copyBytes, buffer, out);
        outBuff[copyBytes + out] = 0;
        currentPos += out + 1;

        return true;
    }

    return false;
}

//----------------------------------------------------------------------
// Reads a line from the file.
//----------------------------------------------------------------------

bool FileReader::getLine(std::string& outString) {
    if (eofReached) {
        return false;
    }

    outString.clear();

    currentPos += leaveSymbols;

    buffer[BUFFER_SIZE] = 0;

    std::size_t out;
    if (getEolPos(buffer + currentPos, fileBytesRead - currentPos, out)) {
        if (out + currentPos >= BUFFER_SIZE) {
            std::cerr << "Wrong file format." << std::endl;
            exit(EXIT_FAILURE);
        }
        outString.append(buffer + currentPos, out);
        currentPos += out + 1;
        return true;
    }
    else {
        // reached end of buffer
        std::size_t copyBytes = fileBytesRead - currentPos;
        outString.append(buffer + currentPos, copyBytes);
        fileBytesRead = fileRead(buffer, sizeof (char), BUFFER_SIZE, file);

        currentPos = 0;

        if (copyBytes == 0 && windowsNewLine) {
            if (buffer[0] == '\n') {
                ++currentPos; // like leaveSymbols == 1
            }
        }

        const bool found = getEolPos(buffer + currentPos, fileBytesRead - currentPos, out);

        if (out + copyBytes >= BUFFER_SIZE - 1) {
            std::cerr << "Wrong file format." << std::endl;
            exit(EXIT_FAILURE);
        }

        if (!found) {
            eofReached = true;

            if (out + copyBytes == 0) {
                return false;
            }
        }

        outString.append(buffer, out);
        currentPos += out + 1;

        return true;
    }

    return false;
}

//----------------------------------------------------------------------
// Stores the input string in the file.
//----------------------------------------------------------------------

bool FileReader::putString(const std::string& buff) {
    if (fileType == RAW) {
        return fwrite(buff.c_str(), sizeof (char), buff.length(), file) == buff.length();
    }
    else if (fileType == GZIP) {
        return static_cast<std::size_t> (gzwrite(gz_file, buff.c_str(), static_cast<unsigned>(sizeof (char) * buff.length()))) == buff.length();
    }
    else {
        std::cerr << "File type not defined." << std::endl;
        exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------
// Returns the current file cursor.
//----------------------------------------------------------------------

long FileReader::tellPos() {
    long filePos = 0;

    if (fileType == RAW) {
        filePos = ftell(file);
    }
    else if (fileType == GZIP) {
        filePos = gztell(gz_file);
    }
    else {
        std::cerr << "File type not defined." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (filePos == 0) {
        return filePos;
    }
    return filePos - static_cast<long>(fileBytesRead) + static_cast<long>(currentPos);
}

//----------------------------------------------------------------------
// Moves the file cursor to the specified position.
//----------------------------------------------------------------------

long FileReader::seekPos(long pos) {
    if (fileType == RAW) {
        return fseek(file, pos, SEEK_SET);
    }
    else if (fileType == GZIP) {
        return gzseek(gz_file, pos, SEEK_SET);
    }
    else {
        std::cerr << "File type not defined." << std::endl;
        exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------
// Checks, whether the EOF have been reached.
//----------------------------------------------------------------------

bool FileReader::eof() const {
    return eofReached;
}

//----------------------------------------------------------------------
// Sets file name and file type (uncompressed, gzipped etc.).
//----------------------------------------------------------------------

void FileReader::setFileName(const std::string& _fileName, FileType _fileType) {
    fileName = _fileName;
    fileType = _fileType;
}

//----------------------------------------------------------------------
// Opens file in the specified mode (for read or for write).
//----------------------------------------------------------------------

bool FileReader::openFile(OpenMode openMode) {
    const char* mode;
    if (openMode == READ) {
        mode = "rb";
    }
    else if (openMode == WRITE) {
        mode = "wb";
    }
    else {
        std::cerr << "Open mode not defined" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (fileType == RAW) {
        file = fopen(fileName.c_str(), mode);
        return file != NULL;
    }
    else if (fileType == GZIP) {
        gz_file = gzopen(fileName.c_str(), mode);
        if (gz_file != NULL) {
            gzbuffer(gz_file, GZIP_BUFFER_SIZE);
        }
        return gz_file != NULL;
    }
    else {
        std::cerr << "File type not defined." << std::endl;
        exit(EXIT_FAILURE);
    }
    return false;
}

//----------------------------------------------------------------------
// Closes the file.
//----------------------------------------------------------------------

void FileReader::close() {
    leaveSymbols = 0;
    currentPos = 0;
    fileBytesRead = 0;
    eofReached = false;
    windowsNewLine = false;
    if (fileType == RAW) {
        if (file != NULL) {
            fclose(file);
            file = NULL;
        }
    }
    else if (fileType == GZIP) {
        if (gz_file != NULL) {
            gzclose(gz_file);
            gz_file = NULL;
        }
    }
}

//----------------------------------------------------------------------
// Performs actual reading.
//----------------------------------------------------------------------

inline std::size_t FileReader::fileRead(char* outBuff, std::size_t elemSize, std::size_t elemCount, FILE* file) {
    if (fileType == RAW) {
        return fread(outBuff, elemSize, elemCount, file);
    }
    else if (fileType == GZIP) {
        return gzread(gz_file, outBuff, static_cast<unsigned>(sizeof (char) * elemCount));
    }
    else {
        std::cerr << "File type not defined." << std::endl;
        exit(EXIT_FAILURE);
    }
}
