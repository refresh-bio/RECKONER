/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
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
        std::cout << "ERROR: Cannot open " << read_file_name << " for read check" << std::endl;
        exit(EXIT_FAILURE);
    }

    // initialize variables
    num_reads = 0;

    max_quality_score = -1000;
    min_quality_score = 1000;

    std::string line;

    // header
#ifdef USE_FILE_READER
    while (f_read.getLine(line)) {
#else
    getline(f_read, line);

    while (!f_read.eof()) {
#endif
        // increment the number of reads
        num_reads++;

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

#ifndef USE_FILE_READER
        // header
        getline(f_read, line);
#endif
    }


    if (min_quality_score < MIN_SCORE) {
        std::cout << "ERROR: Illegal quality score" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (min_quality_score < MIN_64_SCORE) {
        quality_score_offset = PHRED33;
    }
    else {
        quality_score_offset = PHRED64;
    }

    f_read.close();

    create_chunks(c_inst_args, read_file_name, read_file_type, num_reads);

    std::cout << "     Number of reads     : " << num_reads << std::endl;
    std::cout << "     Maximum read length : " << max_read_length << std::endl;
    std::cout << "     Quality score offset: " << quality_score_offset << std::endl;
    std::cout << "     Checking input read files: done" << std::endl << std::endl;

    f_log << "     Number of reads     : " << num_reads << std::endl;
    f_log << "     Maximum read length : " << max_read_length << std::endl;
    f_log << "     Quality score offset: " << quality_score_offset << std::endl;
    f_log << "     Checking input read files: done" << std::endl << std::endl;
}

//----------------------------------------------------------------------
// Reads input file and splits it into chunks for a parallel processing.
//----------------------------------------------------------------------

void C_check_read::create_chunks(const C_arg& c_inst_args, const std::string& read_file_name, const FileReader::FileType read_file_type, std::size_t num_reads) {
    const std::size_t n_chunks = N_CHUNKS_PER_THREAD * c_inst_args.n_threads;
    const std::size_t chunk_size = n_chunks == 1 ? 0 : num_reads / (n_chunks - 1);

#ifdef USE_FILE_READER
    FileReader f_read;
#else
    std::ifstream f_read;
#endif   

    chunks.reserve(n_chunks);
    if (chunk_size == 0) {
#ifdef USE_FILE_READER
        chunks.push_back(0);
#else
        chunks.push_back(f_read.beg);
#endif
    }
    else {
#ifdef USE_FILE_READER
        f_read.setFileName(read_file_name, read_file_type);
        if (!f_read.openFile(FileReader::READ)) {
#else
        f_read.open(read_file_name.c_str());
        if (!f_read.is_open()) {
#endif
            std::cout << "ERROR: Cannot open " << read_file_name << " for chunkify" << std::endl;
            exit(EXIT_FAILURE);
        }

        for (std::size_t it = 0; it < num_reads; ++it) {
            if (it % chunk_size == 0) {
#ifdef USE_FILE_READER
                chunks.push_back(f_read.tellPos());
#else
                chunks.push_back(f_read.tellg());
#endif
            }
            std::string temp;
            for (int i = 0; i < READ_LINES; ++i) {
#ifdef USE_FILE_READER
                f_read.getLine(temp);
#else
                getline(f_read, temp);
#endif
            }
        }
        f_read.close();
    }
}
