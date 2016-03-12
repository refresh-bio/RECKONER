/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#include "PerformTesting.h"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "correct_errors.hpp"
#include "time.hpp"
#include "Log.h"
#include <kmc_api/kmc_file.h>

//----------------------------------------------------------------------
// For every input file calls initial checking.
//----------------------------------------------------------------------

void PerformTasks::perform_check_read(const C_arg& c_inst_args, C_time& c_inst_time, std::vector<C_check_read>& check_inputs) {
    std::ofstream& f_log = Log::get_stream();

    std::cout << "Checking input read files" << std::endl;
    f_log << "Checking input read files" << std::endl;

    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << c_inst_args.read_files_names[it] << std::endl;
        f_log << c_inst_args.read_files_names[it] << std::endl;

        time_t rawtime;

        time(&rawtime);
        c_inst_time.vector_start_check_read_file.push_back(asctime(localtime(&rawtime)));

        C_check_read check_read;
        check_read.check_read_file(c_inst_args, c_inst_args.read_files_names[it], c_inst_args.read_files_types[it]);
        check_inputs.push_back(check_read);

        time(&rawtime);
        c_inst_time.vector_end_check_read_file.push_back(asctime(localtime(&rawtime)));
    }
}

//----------------------------------------------------------------------
// Removes new line symbols from the input string.
//----------------------------------------------------------------------

std::string PerformTasks::remove_new_line(std::string in_string) {
    in_string.erase(std::remove(in_string.begin(), in_string.end(), '\n'), in_string.end());
    return in_string;
}

//----------------------------------------------------------------------
// Prints report and correction times.
//----------------------------------------------------------------------

void PerformTasks::summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time) {
    std::ofstream& f_log = Log::get_stream();
    //--------------------------------------------------
    // stdout
    //--------------------------------------------------
    std::cout << "Running Time" << std::endl;
    std::cout << "     Parsing arguments" << std::endl;
    std::cout << "          Start:" << remove_new_line(c_inst_time.start_parse_args) << std::endl;
    std::cout << "          End  :" << remove_new_line(c_inst_time.end_parse_args) << std::endl;

    std::cout << "     Checking reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << "          " << c_inst_args.read_files_names[it] << std::endl;
        std::cout << "               Start:" << remove_new_line(c_inst_time.vector_start_check_read_file[it]) << std::endl;
        std::cout << "               End  :" << remove_new_line(c_inst_time.vector_end_check_read_file[it]) << std::endl;
    }

    std::cout << "     Correcting errors in reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << "          " << c_inst_args.read_files_names[it] << std::endl;
        std::cout << "               Start:" << remove_new_line(c_inst_time.vector_start_correct_errors_in_reads[it]) << std::endl;
        std::cout << "               End  :" << remove_new_line(c_inst_time.vector_end_correct_errors_in_reads[it]) << std::endl;
    }

    if (c_inst_args.nowrite == false) {
        std::cout << "     Writing corrected reads" << std::endl;
        for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
            std::cout << "          " << c_inst_args.read_files_names[it] << std::endl;
            std::cout << "               Start:" << remove_new_line(c_inst_time.vector_start_write_corrected_reads[it]) << std::endl;
            std::cout << "               End  :" << remove_new_line(c_inst_time.vector_end_write_corrected_reads[it]) << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "The program is successfully completed" << std::endl << std::endl;

    //--------------------------------------------------
    // log file
    //--------------------------------------------------
    f_log << "Running Time" << std::endl;
    f_log << "     Parsing arguments" << std::endl;
    f_log << "          Start:" << remove_new_line(c_inst_time.start_parse_args) << std::endl;
    f_log << "          End  :" << remove_new_line(c_inst_time.end_parse_args) << std::endl;

    f_log << "     Checking reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        f_log << "          " << c_inst_args.read_files_names[it] << std::endl;
        f_log << "               Start:" << remove_new_line(c_inst_time.vector_start_check_read_file[it]) << std::endl;
        f_log << "               End  :" << remove_new_line(c_inst_time.vector_end_check_read_file[it]) << std::endl;
    }

    f_log << "     Correcting errors in reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        f_log << "          " << c_inst_args.read_files_names[it] << std::endl;
        f_log << "               Start:" << remove_new_line(c_inst_time.vector_start_correct_errors_in_reads[it]) << std::endl;
        f_log << "               End  :" << remove_new_line(c_inst_time.vector_end_correct_errors_in_reads[it]) << std::endl;
    }

    if (c_inst_args.nowrite == false) {
        f_log << "     Writing corrected reads" << std::endl;
        for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
            f_log << "          " << c_inst_args.read_files_names[it] << std::endl;
            f_log << "               Start:" << remove_new_line(c_inst_time.vector_start_write_corrected_reads[it]) << std::endl;
            f_log << "               End  :" << remove_new_line(c_inst_time.vector_end_write_corrected_reads[it]) << std::endl;
        }
    }

    f_log << std::endl;
    f_log << "The program is successfully completed" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// Opens KMC database and calls correcting for every input file.
//----------------------------------------------------------------------

void PerformTasks::perform_correction(const C_arg& c_inst_args, C_time& c_inst_time, std::vector<C_check_read>& check_inputs) {
    std::ofstream& f_log = Log::get_stream();

    std::cout << "Correcting errors in reads" << std::endl;
    f_log << "Correcting errors in reads" << std::endl;

    CKMCFile kmc_file;
    if (!kmc_file.OpenForRA(c_inst_args.read_files_names.front() + LIST_FILE_EXTENSION)) {
        std::cerr << "ERROR: cannot open KMC files." << std::endl;
        exit(EXIT_FAILURE);
    }

    uint32 k, mode, counter_size, prefix_length, signature_length, min;
    uint64 max, total;

    kmc_file.Info(k, mode, counter_size, prefix_length, signature_length, min, max, total);

    if (k != c_inst_args.kmer_length) {
        std::cerr << "Error: KMC and provided k-mer lenghts are different." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<C_correct_errors> correct_errors;

    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << c_inst_args.read_files_names[it] << std::endl;
        f_log << c_inst_args.read_files_names[it] << std::endl;

        time_t rawtime;

        // corrent errors in reads
        time(&rawtime);
        c_inst_time.vector_start_correct_errors_in_reads.push_back(asctime(localtime(&rawtime)));
        correct_errors.push_back(C_correct_errors(c_inst_args, check_inputs[it], c_inst_args.read_files_names[it], c_inst_args.read_files_types[it]));
        correct_errors.back().correct_errors_in_reads(c_inst_args, c_inst_args.error_correction_info_files_names[it], kmc_file);

        time(&rawtime);
        c_inst_time.vector_end_correct_errors_in_reads.push_back(asctime(localtime(&rawtime)));
    }
    
    std::cout << "Writing corrected reads into files" << std::endl;
    f_log << "Writing corrected reads into files" << std::endl;
    
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << c_inst_args.read_files_names[it] << std::endl;
        f_log << c_inst_args.read_files_names[it] << std::endl;
        
        time_t rawtime;
        // write corrected reads into output files
        time(&rawtime);
        c_inst_time.vector_start_write_corrected_reads.push_back(asctime(localtime(&rawtime)));
        if (c_inst_args.nowrite == false) {
            correct_errors[it].write_corrected_reads(c_inst_args, c_inst_args.error_correction_info_files_names[it], c_inst_args.corrected_read_files_names[it], c_inst_args.read_files_types[it]);
        }

        time(&rawtime);
        c_inst_time.vector_end_write_corrected_reads.push_back(asctime(localtime(&rawtime)));

        // remove error correction information files
        correct_errors[it].remove_error_correction_info_files(c_inst_args, c_inst_args.error_correction_info_files_names[it]);
    }

    // summarize output results
    summarize_outputs(c_inst_args, c_inst_time);
}
