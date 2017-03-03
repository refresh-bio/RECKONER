/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.0
 * 
 */

#include "RunExternal.h"
#include "PerformTasks.h"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "correct_errors.hpp"
#include "time.hpp"
#include "PrepareKMCDb.h"
#include "DetermineParameters.hpp"
#include <kmc_api/kmc_file.h>



//----------------------------------------------------------------------
// For every input file calls initial checking.
//----------------------------------------------------------------------

void PerformTasks::perform_check_read() {
    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "CHECKING INPUT READ FILES" << std::endl;
    std::cout << "##################################################" << std::endl;
    f_log << std::endl;
    f_log << "##################################################" << std::endl;
    f_log << "CHECKING INPUT READ FILES" << std::endl;
    f_log << "##################################################" << std::endl;

    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << c_inst_args.read_files_names[it] << std::endl;
        f_log << c_inst_args.read_files_names[it] << std::endl;

        Timer check_read_file_time;
        check_read_file_time.startTimer();

        C_check_read check_read;
        check_read.check_read_file(c_inst_args, c_inst_args.read_files_names[it], c_inst_args.read_files_types[it]);
        check_inputs.push_back(check_read);

        check_read_file_time.stopTimer();
        c_inst_time.vector_check_read_file.push_back(check_read_file_time);
    }
}

//----------------------------------------------------------------------
// Determines parameters values (if needed).
//----------------------------------------------------------------------

void PerformTasks::perform_parameter_determination() {
    if (!c_inst_args.is_kmer_length_user_defined) {
        std::cout << std::endl;
        std::cout << "##################################################" << std::endl;
        std::cout << "DETERMINING PARAMETERS" << std::endl;
        std::cout << "##################################################" << std::endl;
        f_log << std::endl;
        f_log << "##################################################" << std::endl;
        f_log << "DETERMINING PARAMETERS" << std::endl;
        f_log << "##################################################" << std::endl;

        Timer determine_parameters_time;
        determine_parameters_time.startTimer();

        DetermineParameters determineParameters(c_inst_args, check_inputs);

        if (!c_inst_args.is_genome_size_user_defined) {
            c_inst_args.genome_size = determineParameters.determineGenomeSize();
            std::cout << "Estimated genome size: " << c_inst_args.genome_size << std::endl;
            f_log << "Estimated genome size: " << c_inst_args.genome_size << std::endl;
        }
        c_inst_args.kmer_length = determineParameters.determineKmerLength();
        std::cout << "Determined k-mer length: " << c_inst_args.kmer_length << std::endl;
        f_log << "Determined k-mer length: " << c_inst_args.kmer_length << std::endl;

        determine_parameters_time.stopTimer();
        c_inst_time.determine_parameters = determine_parameters_time;
    }
}

//----------------------------------------------------------------------
// Calls k-mer counting and KMC database cutting.
//----------------------------------------------------------------------

void PerformTasks::perform_kmcdb_creation() {
    PrepareKMCDb prepareKMCDb(c_inst_args);

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "K - MER COUNTING" << std::endl;
    std::cout << "##################################################" << std::endl;
    f_log << std::endl;
    f_log << "##################################################" << std::endl;
    f_log << "K - MER COUNTING" << std::endl;
    f_log << "##################################################" << std::endl;

    Timer kmer_count_time;
    kmer_count_time.startTimer();

    prepareKMCDb.countKmers();

    kmer_count_time.stopTimer();
    c_inst_time.kmer_count = kmer_count_time;

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "DETERMINING CUTOFF THRESHOLD" << std::endl;
    std::cout << "##################################################" << std::endl;
    f_log << std::endl;
    f_log << "##################################################" << std::endl;
    f_log << "DETERMINING CUTOFF THRESHOLD" << std::endl;
    f_log << "##################################################" << std::endl;

    Timer determine_cutoff_threshold_time;
    determine_cutoff_threshold_time.startTimer();

    unsigned cutoff = prepareKMCDb.determineCutoff();

    determine_cutoff_threshold_time.stopTimer();
    c_inst_time.determine_cutoff_threshold = determine_cutoff_threshold_time;

    std::cout << "Cutoff: " << cutoff << std::endl;
    f_log << "Cutoff: " << cutoff << std::endl;

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "REMOVING UNTRUSTED K - MERS" << std::endl;
    std::cout << "##################################################" << std::endl;
    f_log << std::endl;
    f_log << "##################################################" << std::endl;
    f_log << "REMOVING UNTRUSTED K - MERS" << std::endl;
    f_log << "##################################################" << std::endl;

    Timer remove_untrusted_kmers_time;
    remove_untrusted_kmers_time.startTimer();

    prepareKMCDb.filterKmers(cutoff);

    remove_untrusted_kmers_time.stopTimer();
    c_inst_time.remove_untrusted_kmers = remove_untrusted_kmers_time;

    std::cout << std::endl;

    prepareKMCDb.removeKMCDatabase();
}

//----------------------------------------------------------------------
// Prints report and correction times.
//----------------------------------------------------------------------

void PerformTasks::summarize_outputs() {
    //--------------------------------------------------
    // stdout
    //--------------------------------------------------
    std::cout << "Running Time" << std::endl;

    std::cout << "     Checking reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << "          " << c_inst_args.read_files_names[it] << std::endl;
        std::cout << "               Start:" << c_inst_time.vector_check_read_file[it].getStartTime() << std::endl;
        std::cout << "               End  :" << c_inst_time.vector_check_read_file[it].getStopTime() << std::endl;
    }

    if (c_inst_args.is_kmer_length_user_defined == false) {
        std::cout << "     Determining parameters" << std::endl;
        std::cout << "          Start:" << c_inst_time.determine_parameters.getStartTime() << std::endl;
        std::cout << "          End  :" << c_inst_time.determine_parameters.getStopTime() << std::endl;
    }

    std::cout << "     Counting k-mers" << std::endl;
    std::cout << "          Start:" << c_inst_time.kmer_count.getStartTime() << std::endl;
    std::cout << "          End  :" << c_inst_time.kmer_count.getStopTime() << std::endl;

    std::cout << "     Determining cutoff threshold" << std::endl;
    std::cout << "          Start:" << c_inst_time.determine_cutoff_threshold.getStartTime() << std::endl;
    std::cout << "          End  :" << c_inst_time.determine_cutoff_threshold.getStopTime() << std::endl;

    std::cout << "     Removing untrusted k-mers" << std::endl;
    std::cout << "          Start:" << c_inst_time.remove_untrusted_kmers.getStartTime() << std::endl;
    std::cout << "          End  :" << c_inst_time.remove_untrusted_kmers.getStopTime() << std::endl;

    std::cout << "     Correcting errors in reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        std::cout << "          " << c_inst_args.read_files_names[it] << std::endl;
        std::cout << "               Start:" << c_inst_time.vector_correct_errors_in_reads[it].getStartTime() << std::endl;
        std::cout << "               End  :" << c_inst_time.vector_correct_errors_in_reads[it].getStopTime() << std::endl;
    }

    std::cout << std::endl;
    std::cout << "The program is successfully completed" << std::endl << std::endl;

    //--------------------------------------------------
    // log file
    //--------------------------------------------------
    f_log << "Running Time" << std::endl;

    f_log << "     Checking reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        f_log << "          " << c_inst_args.read_files_names[it] << std::endl;
        f_log << "               Start:" << c_inst_time.vector_check_read_file[it].getStartTime() << std::endl;
        f_log << "               End  :" << c_inst_time.vector_check_read_file[it].getStopTime() << std::endl;
    }

    if (c_inst_args.is_kmer_length_user_defined == false) {
        f_log << "     Determining parameters" << std::endl;
        f_log << "          Start:" << c_inst_time.determine_parameters.getStartTime() << std::endl;
        f_log << "          End  :" << c_inst_time.determine_parameters.getStopTime() << std::endl;
    }

    f_log << "     Counting k-mers" << std::endl;
    f_log << "          Start:" << c_inst_time.kmer_count.getStartTime() << std::endl;
    f_log << "          End  :" << c_inst_time.kmer_count.getStopTime() << std::endl;

    f_log << "     Determining cutoff threshold" << std::endl;
    f_log << "          Start:" << c_inst_time.determine_cutoff_threshold.getStartTime() << std::endl;
    f_log << "          End  :" << c_inst_time.determine_cutoff_threshold.getStopTime() << std::endl;

    f_log << "     Removing untrusted k-mers" << std::endl;
    f_log << "          Start:" << c_inst_time.remove_untrusted_kmers.getStartTime() << std::endl;
    f_log << "          End  :" << c_inst_time.remove_untrusted_kmers.getStopTime() << std::endl;

    f_log << "     Correcting errors in reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        f_log << "          " << c_inst_args.read_files_names[it] << std::endl;
        f_log << "               Start:" << c_inst_time.vector_correct_errors_in_reads[it].getStartTime() << std::endl;
        f_log << "               End  :" << c_inst_time.vector_correct_errors_in_reads[it].getStopTime() << std::endl;
    }

    f_log << std::endl;
    f_log << "The program is successfully completed" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// Opens KMC database and calls correcting for every input file.
//----------------------------------------------------------------------

void PerformTasks::perform_correction() {
    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "CORRECTING ERRORS" << std::endl;
    std::cout << "##################################################" << std::endl;
    Log::get_stream() << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    Log::get_stream() << "CORRECTING ERRORS" << std::endl;
    Log::get_stream() << "##################################################" << std::endl;

    CKMCFile kmc_file;
    if (!kmc_file.OpenForRA(c_inst_args.kmc_filtered_database_name)) {
        std::cerr << "ERROR: cannot open KMC files." << std::endl;
        exit(EXIT_FAILURE);
    }

    uint32 k, mode, counter_size, prefix_length, signature_length, min;
    uint64 max, total;

    kmc_file.Info(k, mode, counter_size, prefix_length, signature_length, min, max, total);

    if (k != c_inst_args.kmer_length) {
        std::cerr << "Error: KMC and provided k-mer lengths are different." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (std::size_t it = 0; it < c_inst_args.read_files_names.size(); ++it) {
        C_correct_errors reads_corrector(c_inst_args, c_inst_args.error_correction_info_files_names[it], c_inst_args.read_files_names[it], c_inst_args.corrected_read_files_names[it], c_inst_args.read_files_types[it], check_inputs[it]);

        std::cout << c_inst_args.read_files_names[it] << std::endl;
        f_log << c_inst_args.read_files_names[it] << std::endl;

        Timer correct_errors_in_reads_timer;
        correct_errors_in_reads_timer.startTimer();

        // corrent errors in reads
        reads_corrector.correct_errors_in_reads(kmc_file);

        correct_errors_in_reads_timer.stopTimer();
        c_inst_time.vector_correct_errors_in_reads.push_back(correct_errors_in_reads_timer);
    }
}

//----------------------------------------------------------------------
// Prints summary and removes remaining KMC files.
//----------------------------------------------------------------------

void PerformTasks::finalize() {
    // summarize output results
    summarize_outputs();

    PrepareKMCDb prepareKMCDb(c_inst_args);
    prepareKMCDb.removeFilteredKMCDatabase();
}