/*
 * RECKONER - Read Error Corrector Based on KMC
 *
 * This software is distributed under GNU GPL 3 license.
 *
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.0
 *
*/

#ifndef _FASTQ_READER_H
#define _FASTQ_READER_H

#include "FileReader.h"
#include "Log.h"
#include <zlib.h>
#include <string>
#include <stack>
#include <queue>
#include <mutex>
#include <condition_variable>



class MemoryPool {
    std::mutex bufferAvailableMutex;
    std::condition_variable buffersAvailable;
    std::stack<char*> pool;

public:
    MemoryPool(std::size_t bufferSize, std::size_t poolSize);
    ~MemoryPool();

    char* getBuffer();
    void releaseBuffer(char* buffer);
};



//----------------------------------------------------------------------
// CFastqReader - based on KMC 2.3.0
// Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot
//----------------------------------------------------------------------

class FastqReader {
    std::ofstream& f_log;

    MemoryPool *pmm_fastq;

    std::string input_file_name;
    FileReader::FileType mode;

    FILE *in;
    gzFile_s *in_gzip;
    int bzerror;

    std::size_t part_size;

    char *part;
    std::size_t part_filled;

    unsigned gzip_buffer_size;

    bool SkipNextEOL(char *part, long long &pos, long long max_pos);

    bool IsEof();

public:
    FastqReader(MemoryPool *_pmm_fastq, unsigned _gzip_buffer_size);

    ~FastqReader();

    static std::size_t OVERHEAD_SIZE;

    bool SetNames(std::string _input_file_name, FileReader::FileType _file_type);

    bool SetPartSize(std::size_t _part_size);

    bool OpenFiles();

    bool GetPart(char *&_part, std::size_t &_size);
};



class ReadsChunk {
private:
    std::size_t pos;
    MemoryPool* memoryPool;

    char* buffer;
    std::size_t size;

    void setChunk(char* _buffer, std::size_t bufferSize, MemoryPool* _memoryPool);
    void releaseMemory();

public:
    ReadsChunk() : pos(0), memoryPool(NULL), buffer(NULL), size(0) {}

    ~ReadsChunk() {
        releaseMemory();
    }
    bool getLine(std::string& line);

    friend class FastqReaderWrapper;
};



class FastqReaderWrapper {
    MemoryPool memoryPool;
    FastqReader reader;

    bool finished;

    std::mutex partsAvailableMutex;
    std::condition_variable partsAvailable;
    std::queue<std::pair<char*, size_t>> partsQueue;

public:
    FastqReaderWrapper(const FastqReaderWrapper&) = delete;
    FastqReaderWrapper(std::size_t bufferSize, std::size_t poolSize);

    bool openFile(std::string inputFileName, FileReader::FileType fileType);

    void operator()();

    bool getChunk(ReadsChunk& readsChunk);
};

#endif
