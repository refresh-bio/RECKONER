/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.2
 * 
 */

#include "check_inputs.hpp"
#include "FileReader.h"
#include "Log.h"
#include <string>



//----------------------------------------------------------------------
// Reads and analyzes entire input file.
//----------------------------------------------------------------------

void C_check_read::check_read_file(const C_arg& c_inst_args, const std::string& read_file_name, const FileReader::FileType read_file_type) {
    std::ofstream& f_log = Log::get_stream();
    // open input read files
#ifdef USE_FILE_READER
    FileReader f_read;
    f_read.setFileName(read_file_name, read_file_type);
    if (!f_read.openFile(FileReader::READ)) {
#else
    std::ifstream f_read(read_file_name.c_str());
    if (!f_read.is_open()) {
#endif
        std::cerr << "ERROR: Cannot open " << read_file_name << " for read check" << std::endl;
        exit(EXIT_FAILURE);
    }

    // initialize variables
    num_reads = 0;

    max_quality_score = -1000;
    min_quality_score = 1000;

    std::string line;

#ifdef USE_FILE_READER
    chunks.push_back(f_read.tellPos());
#else
    chunks.push_back(f_read.tellg());
#endif

    // header
#ifdef USE_FILE_READER
    while (f_read.getLine(line)) {
#else
    getline(f_read, line);

    while (!f_read.eof()) {
#endif

        // sequence
#ifdef USE_FILE_READER
        f_read.getLine(line);
#else
        getline(f_read, line);
#endif

        std::size_t read_length = line.length();

        if (read_length > max_read_length) {
            max_read_length = read_length;
        }

        // connector
#ifdef USE_FILE_READER
        f_read.getLine(line);
#else
        getline(f_read, line);
#endif

        // quality score
#ifdef USE_FILE_READER
        f_read.getLine(line);
#else
        getline(f_read, line);
#endif
        for (std::size_t it = 0; it < line.length(); it++) {
            if ((int) line[it] > max_quality_score) {
                max_quality_score = (int) line[it];
            }

            if ((int) line[it] < min_quality_score) {
                min_quality_score = (int) line[it];
            }
        }

        // increment the number of reads
        num_reads++;

        if (num_reads % chunk_size == 0) {
#ifdef USE_FILE_READER
            chunks.push_back(f_read.tellPos());
#else
            chunks.push_back(f_read.tellg());
#endif
        }

#ifndef USE_FILE_READER
        // header
        getline(f_read, line);
#endif
    }

    last_chunk_size = num_reads % chunk_size;

    // no elements in last chunk, so chunks.back() indicates end of the file
    if (last_chunk_size == 0) {
        chunks.pop_back();

        // now the last chunk is full, not empty
        last_chunk_size = chunk_size;
    } else if (num_reads == 0) {
        chunks.clear();
    }

    f_read.close();

    if (min_quality_score < MIN_SCORE) {
        std::cerr << "ERROR: Illegal quality score" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    } else if (min_quality_score < MIN_64_SCORE) {
        quality_score_offset = PHRED33;
    } else {
        quality_score_offset = PHRED64;
    }

    std::cout << "     Number of reads     : " << num_reads << std::endl;
    std::cout << "     Maximum read length : " << max_read_length << std::endl;
    std::cout << "     Quality score offset: " << quality_score_offset << std::endl;
    std::cout << "     Checking input read files: done" << std::endl;

    f_log << "     Number of reads     : " << num_reads << std::endl;
    f_log << "     Maximum read length : " << max_read_length << std::endl;
    f_log << "     Quality score offset: " << quality_score_offset << std::endl;
    f_log << "     Checking input read files: done" << std::endl;
}
