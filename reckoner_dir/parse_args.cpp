/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#include "parse_args.hpp"
#include "Log.h"
#include <string.h>



//----------------------------------------------------------------------
// Parses program arguments.
//----------------------------------------------------------------------

void C_arg::read_args() {
    std::cout << "RECKONER ver. " << VERSION << std::endl;
    if (num_args == 1) {
        print_usage();
        exit(EXIT_SUCCESS);
    }

    //--------------------------------------------------
    // parse arguments
    //--------------------------------------------------
    for (int it_arg = 1; it_arg < num_args; it_arg++) {
        if (strcmp(args[it_arg], "-help") == 0) {
            print_usage();
            exit(EXIT_SUCCESS);
        }
        else if (strcmp(args[it_arg], "-read") == 0) {
            if (it_arg <= num_args - 2) {
                read_files_names.push_back(args[it_arg + 1]);
                it_arg++;
            }
            else {
                std::cout << std::endl << "ERROR: The read file name is not specified" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-prefix") == 0) {
            if (it_arg <= num_args - 2) {
                prefix = args[it_arg + 1];
                it_arg++;
            }
            else {
                std::cout << std::endl << "ERROR: The prefix is not specified" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-kmerlength") == 0) {
            if (it_arg <= num_args - 2) {
                kmer_length = atoi(args[it_arg + 1]);
                it_arg++;
            }
            else {
                std::cout << std::endl << "ERROR: The k-mer length is not specified" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-extend") == 0) {
            if (it_arg <= num_args - 2) {
                extend = atoi(args[it_arg + 1]);
                it_arg++;
            }
            else {
                std::cout << std::endl << "ERROR: The max extension is not specified" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-threads") == 0) {
            if (it_arg <= num_args - 2) {
                n_threads = atoi(args[it_arg + 1]);
                it_arg++;
            }
            else {
                std::cout << std::endl << "ERROR: The number of threads is not specified" << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-nowrite") == 0) {
            nowrite = true;
        }
        else {
            std::cout << std::endl << "ERROR: Illegal option " << args[it_arg] << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }
    }


    //--------------------------------------------------
    // check options
    //--------------------------------------------------
    if (read_files_names.empty()) {
        std::cout << std::endl << "ERROR: No read file name is specified." << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (prefix.empty()) {
        prefix = ".";
    }

    if (kmer_length == 0) {
        std::cout << std::endl << "ERROR: k-mer length not specified" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (kmer_length < MIN_KMER_LENGTH) {
        std::cout << std::endl << "ERROR: k-mer length should be >= " << MIN_KMER_LENGTH << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (extend < 1) {
        std::cout << std::endl << "ERROR: The max extension should be >= 1" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    //--------------------------------------------------
    // file names
    //--------------------------------------------------

    read_files_types.resize(read_files_names.size(), FileReader::NOT_DEFINED);

    for (std::size_t it = 0; it < read_files_names.size(); ++it) {
        std::ifstream f_tmp;
        f_tmp.open(read_files_names[it].c_str());

        if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << read_files_names[it] << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }

        f_tmp.close();

        std::size_t slash_pos = read_files_names[it].find_last_of("/");
        std::string input_file_name(slash_pos == std::string::npos ? read_files_names[it] : read_files_names[it].substr(slash_pos + 1));
        std::string last_extension, pre_last_extension;

        std::size_t dot_pos = input_file_name.find_last_of(".");
        if (dot_pos != std::string::npos) {
            last_extension = input_file_name.substr(dot_pos);
            input_file_name = input_file_name.substr(0, dot_pos);

            dot_pos = input_file_name.find_last_of(".");
            if (dot_pos != std::string::npos) {
                pre_last_extension = input_file_name.substr(dot_pos);
                input_file_name = input_file_name.substr(0, dot_pos);
            }
        }

        if (log_file_name.empty()) {
            // set a log file name
            log_file_name = prefix + "/" + input_file_name + LOG_FILE_EXTENSION;
        }

        // set error_correction information file names
        error_correction_info_files_names.push_back(prefix + "/" + input_file_name + CORRECTION_FILE_EXTENSION);

        // set error_corrected read file names
        corrected_read_files_names.push_back(prefix + "/" + input_file_name + CORRECTED_FILE_EXTENSION + pre_last_extension + last_extension);

        std::transform(last_extension.begin(), last_extension.end(), last_extension.begin(), ::tolower);
        if (last_extension == GZ_FILE_EXTENSION) {
            read_files_types[it] = FileReader::GZIP;
        }
        else {
            read_files_types[it] = FileReader::RAW;
        }
    }

    Log::open_log_file(log_file_name);

    //--------------------------------------------------
    // print options
    //--------------------------------------------------
    std::ofstream& f_log = Log::get_stream();

    std::cout << std::endl;
    std::cout << "Parsing arguments is finished" << std::endl;

    for (std::size_t it = 0; it < read_files_names.size(); ++it) {
        std::cout << "     Read File " << (it + 1) << std::endl;
        std::cout << "          Name:                    : " << read_files_names[it] << std::endl;
        std::string compressed = "N";
        if (read_files_types[it] == FileReader::GZIP) {
            compressed = "Y (GZIP)";
        }
        std::cout << "          Compressed:              : " << compressed << std::endl;
    }
    std::cout << "     Log File Name                 : " << log_file_name << std::endl;
    std::cout << "     K-mer Length                  : " << kmer_length << std::endl;
    std::cout << "     Number of Threads             : " << n_threads << std::endl;
    std::cout << std::endl;

    f_log << "RECKONER ver. " << VERSION << std::endl;
    f_log << std::endl;
    f_log << "Parsing arguments is finished" << std::endl;
    for (std::size_t it = 0; it < read_files_names.size(); ++it) {
        f_log << "     Read File " << (it + 1) << std::endl;
        f_log << "          Name:                    : " << read_files_names[it] << std::endl;
        std::string compressed = "N";
        if (read_files_types[it] == FileReader::GZIP) {
            compressed = "Y (GZIP)";
        }
        f_log << "          Compressed:              : " << compressed << std::endl;
    }
    f_log << "     Log File Name                 : " << log_file_name << std::endl;
    f_log << "     K-mer Length                  : " << kmer_length << std::endl;
    f_log << "     Number of Threads             : " << n_threads << std::endl;
    f_log << std::endl;
}

//----------------------------------------------------------------------
// Prints arguments usage.
//----------------------------------------------------------------------

void C_arg::print_usage() {
    std::cout << std::endl;
    std::cout << "USAGE: " << args[0] << " <ARGUMENTS>" << std::endl;
    std::cout << std::endl;
    std::cout << "-help                   - print this text" << std::endl;
    std::cout << "-read file              - FASTQ read file name, can be passed many times" << std::endl;
    std::cout << "-prefix                 - output directory, default current directory (.)" << std::endl;
    std::cout << "-kmerlength             - length of k-mers" << std::endl;
    std::cout << "-extend                 - max extend length, default 2" << std::endl;
    std::cout << "-threads                - number of correcting threads, default number of available virtual cores" << std::endl;
    std::cout << "-nowrite                - no output read" << std::endl;

    std::cout << std::endl;
}
