/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 2.1
*
*/

#include <thread>
#include "FastqReader.h"
#include "define.hpp"
#include <cstring>
#include <algorithm>
#include <iostream>



std::size_t FastqReader::OVERHEAD_SIZE = 1 << 16;



//----------------------------------------------------------------------
// Allocates memory for parts.
//----------------------------------------------------------------------

MemoryPool::MemoryPool(std::size_t bufferSize, std::size_t poolSize) {
    for (std::size_t i = 0; i < poolSize; i++) {
        pool.push(new char[bufferSize]);
    }
}



//----------------------------------------------------------------------
// Realeases parts memory.
//----------------------------------------------------------------------

MemoryPool::~MemoryPool() {
    const std::size_t poolSize = pool.size();
    for (std::size_t i = 0; i < poolSize; i++) {
        delete[] pool.top();
        pool.pop();
    }
}



//----------------------------------------------------------------------
// Waits for memory availability and returns one.
//----------------------------------------------------------------------

char* MemoryPool::getBuffer() {
    std::unique_lock<std::mutex> lck(bufferAvailableMutex);

    buffersAvailable.wait(lck, [this] { return !pool.empty(); });

    char* result = pool.top();
    pool.pop();
    return result;
}



//----------------------------------------------------------------------
// Releases memory to pool.
//----------------------------------------------------------------------

void MemoryPool::releaseBuffer(char* buffer) {
    std::unique_lock<std::mutex> lck(bufferAvailableMutex);

    pool.push(buffer);
    buffersAvailable.notify_one();
}



bool FastqReader::SkipNextEOL(char *part, long long &pos, long long max_pos) {
    long long i;
    for (i = pos; i < max_pos - 2; ++i)
        if ((part[i] == '\n' || part[i] == '\r') && !(part[i + 1] == '\n' || part[i + 1] == '\r'))
            break;

    if (i >= max_pos - 2)
        return false;

    pos = i + 1;

    return true;
}



bool FastqReader::IsEof() {
    if (mode == FileReader::RAW)
        return feof(in) != 0;
    else if (mode == FileReader::GZIP)
        return gzeof(in_gzip) != 0;

    return true;
}



FastqReader::FastqReader(MemoryPool *_pmm_fastq, unsigned _gzip_buffer_size) : c_err(std::cerr) {
    pmm_fastq = _pmm_fastq;

    // Input file mode (default: uncompressed)
    mode = FileReader::RAW;

    // Pointers to input files in various formats (uncompressed, gzip-compressed, bzip2-compressed)
    in = NULL;
    in_gzip = NULL;

    // Size and pointer for the buffer
    part_size = 1 << 23;
    part = NULL;

    gzip_buffer_size = _gzip_buffer_size;
}



FastqReader::~FastqReader() {
    if (mode == FileReader::RAW)
    {
        if (in)
            fclose(in);
    }
    else if (mode == FileReader::GZIP)
    {
        if (in_gzip)
            gzclose(in_gzip);
    }

    if (part)
        pmm_fastq->releaseBuffer(part);
}



bool FastqReader::SetPartSize(std::size_t _part_size) {
    if (in || in_gzip)
        return false;

    if (_part_size < (1 << 20) || _part_size >(1 << 30))
        return false;

    part_size = _part_size;

    return true;
}


#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4996)
#endif
bool FastqReader::OpenFiles() {
    if (in || in_gzip)
        return false;

    // Uncompressed file
    if (mode == FileReader::RAW)
    {
        if ((in = fopen(input_file_name.c_str(), "rb")) == NULL)
            return false;
    }
    // Gzip-compressed file
    else if (mode == FileReader::GZIP)
    {
        if ((in_gzip = gzopen(input_file_name.c_str(), "rb")) == NULL)
            return false;
        gzbuffer(in_gzip, gzip_buffer_size);
    }

    // Reserve via PMM
    part = pmm_fastq->getBuffer();

    part_filled = 0;

    return true;
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif



bool FastqReader::GetPart(char *&_part, std::size_t &_size) {
    if (!in && !in_gzip)
        return false;


    if (IsEof())
        return false;
    std::size_t readed;

    // Read data
    if (mode == FileReader::RAW)
        readed = fread(part + part_filled, 1, part_size - part_filled, in);
    else if (mode == FileReader::GZIP)
        readed = gzread(in_gzip, part + part_filled, (int)(part_size - part_filled));
    else
        readed = 0;				// Never should be here

    long long total_filled = part_filled + readed;
    long long i;

    if (part_filled >= OVERHEAD_SIZE)
    {
        c_err << std::endl << "ERROR: Wrong input file" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (IsEof())
    {
        _part = part;
        _size = total_filled;

        part = NULL;
        return true;
    }

    // Looking for a FASTQ record at the end of the area
    long long line_start[9];
    int j;

    i = total_filled - OVERHEAD_SIZE / 2;
    for (j = 0; j < 9; ++j)
    {
        if (!SkipNextEOL(part, i, total_filled))
            break;
        line_start[j] = i;
    }

    _part = part;
    if (j < 9)
        _size = 0;
    else
    {
        int k;
        for (k = 0; k < 4; ++k)
        {
            if (part[line_start[k] + 0] == '@' && part[line_start[k + 2] + 0] == '+')
            {
                if (part[line_start[k + 2] + 1] == '\n' || part[line_start[k + 2] + 1] == '\r')
                    break;
                if (line_start[k + 1] - line_start[k] == line_start[k + 3] - line_start[k + 2] &&
                    memcmp(part + line_start[k] + 1, part + line_start[k + 2] + 1, line_start[k + 3] - line_start[k + 2] - 1) == 0)
                    break;
            }
        }

        if (k == 4)
            _size = 0;
        else
            _size = line_start[k];
    }

    // Allocate new memory for the buffer

    part = pmm_fastq->getBuffer();
    std::copy(_part + _size, _part + total_filled, part);
    part_filled = total_filled - _size;

    return true;
}



//----------------------------------------------------------------------
// Sets file name.
//----------------------------------------------------------------------

bool FastqReader::SetNames(std::string _input_file_name, FileReader::FileType _file_type) {
    input_file_name = _input_file_name;
    mode = _file_type;

    return true;
}



//----------------------------------------------------------------------
// Set passed memory.
//----------------------------------------------------------------------

void ReadsChunk::setChunk(char* _buffer, std::size_t bufferSize, std::size_t _chunkNo, MemoryPool* _memoryPool) {
    if (buffer != NULL) {
        memoryPool->releaseBuffer(buffer);
    }
    buffer = _buffer;
    size = bufferSize;
    chunkNo = _chunkNo;
    pos = 0;

    memoryPool = _memoryPool;
}



//----------------------------------------------------------------------
// Returns memory to the pool.
//----------------------------------------------------------------------

void ReadsChunk::releaseMemory() {
    if (buffer != NULL) {
        memoryPool->releaseBuffer(buffer);
    }
    buffer = NULL;
}



//----------------------------------------------------------------------
// Reads a single line from chunk.
//----------------------------------------------------------------------

bool ReadsChunk::getLine(std::string& line) {
    if (pos >= size || memoryPool == NULL) {
        return false;
    }

    std::size_t i = pos;
    for (; i < size; ++i) {
        if (buffer[i] == '\n' || buffer[i] == '\r') {
            break;
        }
    }

    line.clear();
    line.append(buffer + pos, i - pos);

    pos = i + 1;
    if (pos < size - 1 && (buffer[pos] == '\n' || buffer[pos] == '\r')) {
        ++pos;
    }

    return true;
}



//----------------------------------------------------------------------
// Initializes the reader.
//----------------------------------------------------------------------

FastqReaderWrapper::FastqReaderWrapper(std::size_t bufferSize, std::size_t poolSize) :
    memoryPool(PART_SIZE, poolSize),       //
    reader(&memoryPool, GZIP_BUFFER_SIZE), // Destruction order is important.
    finished(true) {
    reader.SetPartSize(PART_SIZE);
}



//----------------------------------------------------------------------
// Calls fastq reader file opening.
//----------------------------------------------------------------------

bool FastqReaderWrapper::openFile(std::string _input_file_name, FileReader::FileType _fileType) {
    reader.SetNames(_input_file_name, _fileType);

    if (reader.OpenFiles()) {
        finished = false;
        return true;
    }
    return false;
}



//----------------------------------------------------------------------
// Body of the fastq reader thread.
//----------------------------------------------------------------------

void FastqReaderWrapper::operator()() {
    if (finished) {
        return;
    }

    char* buffer;
    size_t bufferSize;
    size_t chunkNo = 0;
    while (reader.GetPart(buffer, bufferSize)) {
        std::unique_lock<std::mutex> lck(partsAvailableMutex);
        partsQueue.push(PartIndicator(buffer, bufferSize, chunkNo));
        partsAvailable.notify_one();
        ++chunkNo;
    }
    finished = true;
    partsAvailable.notify_all();
}



//----------------------------------------------------------------------
// Waits for parts availability and creates chunk.
//----------------------------------------------------------------------

bool FastqReaderWrapper::getChunk(ReadsChunk& readsChunk) {
    char* buffer;
    size_t bufferSize;
    size_t chunkNo;

    readsChunk.releaseMemory();

    { // for unique_lock
        std::unique_lock<std::mutex> lck(partsAvailableMutex);
        partsAvailable.wait(lck, [this] { return !partsQueue.empty() || finished; });
        if (!partsQueue.empty()) {
            std::tie(buffer, bufferSize, chunkNo) = partsQueue.front();
            partsQueue.pop();
        }
        else {
            return false;
        }
    }

    readsChunk.setChunk(buffer, bufferSize, chunkNo, &memoryPool);

    return true;
}
