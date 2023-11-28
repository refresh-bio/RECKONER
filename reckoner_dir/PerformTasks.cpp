/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.1
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
#include <algorithm>



//----------------------------------------------------------------------
// For every input file calls initial checking.
//----------------------------------------------------------------------

void PerformTasks::perform_check_read() {
    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "CHECKING READS" << std::endl;
    c_log << "##################################################" << std::endl;

    for (std::size_t it = 0; it < c_inst_args.read_files_data.size(); ++it) {
        c_log << c_inst_args.read_files_data[it].input_name << std::endl;

        Timer check_read_file_time;
        check_read_file_time.startTimer();

        C_check_read check_read;
        check_read.check_read_file(c_inst_args.read_files_data[it]);
        check_inputs.push_back(check_read);

        check_read_file_time.stopTimer();
        c_inst_time.vector_check_read_file.push_back(check_read_file_time);
    }
}

//----------------------------------------------------------------------
// Determines parameters values (if needed).
//----------------------------------------------------------------------

void PerformTasks::perform_parameter_determination() {
    if (c_inst_args.is_kmer_length_user_defined) {
        return;
    }

    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "DETERMINING PARAMETERS" << std::endl;
    c_log << "##################################################" << std::endl;

    Timer determine_parameters_time;
    determine_parameters_time.startTimer();

    DetermineParameters determineParameters(c_inst_args, check_inputs);

    if (!c_inst_args.is_genome_size_user_defined) {
        c_inst_args.genome_size = determineParameters.determineGenomeSize();
        c_log << "Estimated genome size: " << c_inst_args.genome_size << std::endl;
    }
    
    c_inst_args.kmer_length = determineParameters.determineKmerLength();
    if (c_inst_args.use_long_kmer) {
        c_inst_args.long_kmer_length = static_cast<std::size_t>(c_inst_args.long_kmer_ratio * c_inst_args.kmer_length);
    }
    else {
        c_inst_args.long_kmer_length = 0;
    }

    c_inst_args.kmer_length = determineParameters.determineKmerLength();
    c_log << "Determined k-mer length: " << c_inst_args.kmer_length << std::endl;

    if (c_inst_args.use_long_kmer) {
        c_log << "Determined long k-mer length: " << c_inst_args.long_kmer_length << std::endl;
    }

    determine_parameters_time.stopTimer();
    c_inst_time.determine_parameters = determine_parameters_time;
}

//----------------------------------------------------------------------
// Calls k-mer counting and KMC database cutting.
//----------------------------------------------------------------------

void PerformTasks::perform_kmcdb_creation() {
    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "K - MER COUNTING" << std::endl;
    c_log << "##################################################" << std::endl;



    if (c_inst_args.reuse_kmc_db) {
        bool db_exists = false;
        bool db_to_reuse = false;

        CKMCFile kmc_file;
        uint32 k, mode, counter_size, prefix_length, signature_length, min;
        uint64 max, total;

        if (kmc_file.OpenForListing(c_inst_args.kmc_filtered_database_name)) {
            db_exists = true;

            kmc_file.Info(k, mode, counter_size, prefix_length, signature_length, min, max, total);
            if (k == c_inst_args.kmer_length) {
                db_to_reuse = true;
            }
        }
        kmc_file.Close();

        if (c_inst_args.use_long_kmer) {
            bool long_db_exists = false;
            bool long_db_to_reuse = false;

            if (kmc_file.OpenForListing(c_inst_args.kmc_long_database_name)) {
                long_db_exists = true;

                kmc_file.Info(k, mode, counter_size, prefix_length, signature_length, min, max, total);
                if (k == c_inst_args.long_kmer_length) {
                    long_db_to_reuse = true;
                }
            }

            if (db_exists && long_db_exists) {
                if (db_to_reuse && long_db_to_reuse) {
                    c_log << "Both KMC databases exist and can be reused." << std::endl;
                    return;
                }
                else {
                    c_log << "At least one of the KMC databases exists, but has different k-mer length. Both of them will be recreated." << std::endl << std::endl;
                }
            }
            else {
                c_log << "At least one of the KMC databases does not exist. Both of them will be created." << std::endl << std::endl;
            }
        }
        else {
            if (db_exists) {
                if (db_to_reuse) {
                    c_log << "KMC database exists and can be reused." << std::endl;
                    return;
                }
                else {
                    c_log << "KMC database exists, but has different k-mer length. It will be recreated." << std::endl << std::endl;
                }
            }
            else {
                c_log << "KMC database does not exist. It will be created." << std::endl << std::endl;
            }
        }
    }

    PrepareKMCDb prepareKMCDb(c_inst_args);

    Timer kmer_count_time;
    kmer_count_time.startTimer();

    prepareKMCDb.countKmers();

    kmer_count_time.stopTimer();
    c_inst_time.kmer_count = kmer_count_time;

    if (c_inst_args.use_long_kmer) {
        Timer long_kmer_count_time;
        long_kmer_count_time.startTimer();

        prepareKMCDb.countLongKmers();

        long_kmer_count_time.stopTimer();
        c_inst_time.long_kmer_count = long_kmer_count_time;
    }



    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "DETERMINING CUTOFF THRESHOLD" << std::endl;
    c_log << "##################################################" << std::endl;

    Timer determine_cutoff_threshold_time;
    determine_cutoff_threshold_time.startTimer();

    unsigned cutoff = prepareKMCDb.determineCutoff();

    determine_cutoff_threshold_time.stopTimer();
    c_inst_time.determine_cutoff_threshold = determine_cutoff_threshold_time;

    c_log << "Cutoff: " << cutoff << std::endl;

    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "REMOVING UNTRUSTED K - MERS" << std::endl;
    c_log << "##################################################" << std::endl;

    Timer remove_untrusted_kmers_time;
    remove_untrusted_kmers_time.startTimer();

    prepareKMCDb.filterKmers(cutoff);

    remove_untrusted_kmers_time.stopTimer();
    c_inst_time.remove_untrusted_kmers = remove_untrusted_kmers_time;

    prepareKMCDb.removeKMCDatabase();
}

//----------------------------------------------------------------------
// Prints report and correction times.
//----------------------------------------------------------------------

void PerformTasks::summarize_outputs() {
    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "RUNNING TIME" << std::endl;
    c_log << "##################################################" << std::endl;

    c_log << "Checking reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_data.size(); ++it) {
        c_log << "     " << c_inst_args.read_files_data[it].input_name << std::endl;
        c_log << "          Start: " << c_inst_time.vector_check_read_file[it].getStartTime() << std::endl;
        c_log << "          End  : " << c_inst_time.vector_check_read_file[it].getStopTime() << std::endl;
    }

    if (c_inst_args.is_kmer_length_user_defined == false) {
        c_log << "Determining parameters" << std::endl;
        c_log << "     Start: " << c_inst_time.determine_parameters.getStartTime() << std::endl;
        c_log << "     End  : " << c_inst_time.determine_parameters.getStopTime() << std::endl;
    }

    if (c_inst_time.kmer_count.wasTimeMeasured()) {
        c_log << "Counting k-mers" << std::endl;
        c_log << "     Start: " << c_inst_time.kmer_count.getStartTime() << std::endl;
        c_log << "     End  : " << c_inst_time.kmer_count.getStopTime() << std::endl;
    }

    if (c_inst_time.long_kmer_count.wasTimeMeasured()) {
        c_log << "Counting long k-mers" << std::endl;
        c_log << "     Start: " << c_inst_time.kmer_count.getStartTime() << std::endl;
        c_log << "     End  : " << c_inst_time.kmer_count.getStopTime() << std::endl;
    }

    if (c_inst_time.determine_cutoff_threshold.wasTimeMeasured()) {
        c_log << "Determining cutoff threshold" << std::endl;
        c_log << "     Start: " << c_inst_time.determine_cutoff_threshold.getStartTime() << std::endl;
        c_log << "     End  : " << c_inst_time.determine_cutoff_threshold.getStopTime() << std::endl;
    }

    if (c_inst_time.remove_untrusted_kmers.wasTimeMeasured()) {
        c_log << "Removing untrusted k-mers" << std::endl;
        c_log << "     Start: " << c_inst_time.remove_untrusted_kmers.getStartTime() << std::endl;
        c_log << "     End  : " << c_inst_time.remove_untrusted_kmers.getStopTime() << std::endl;
    }

    c_log << "Correcting errors in reads" << std::endl;
    for (std::size_t it = 0; it < c_inst_args.read_files_data.size(); ++it) {
        c_log << "     " << c_inst_args.read_files_data[it].input_name << std::endl;
        c_log << "          Start: " << c_inst_time.vector_correct_errors_in_reads[it].getStartTime() << std::endl;
        c_log << "          End  : " << c_inst_time.vector_correct_errors_in_reads[it].getStopTime() << std::endl;
    }

    c_log << std::endl;
    c_log << "The program is successfully completed" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// Opens KMC database and calls correcting for every input file.
//----------------------------------------------------------------------

void PerformTasks::perform_correction() {
    c_log << std::endl;
    c_log << "##################################################" << std::endl;
    c_log << "CORRECTING ERRORS" << std::endl;
    c_log << "##################################################" << std::endl;

    CKMCFile kmc_file, kmc_long_file;
    if (!kmc_file.OpenForRA(c_inst_args.kmc_filtered_database_name)) {
        c_err << "ERROR: cannot open KMC files." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (c_inst_args.use_long_kmer) {
        if (!kmc_long_file.OpenForRA(c_inst_args.kmc_long_database_name)) {
            c_err << "ERROR: cannot open long k-mers KMC files." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    uint32 k, mode, counter_size, prefix_length, signature_length, min;
    uint64 max, total;

    kmc_file.Info(k, mode, counter_size, prefix_length, signature_length, min, max, total);

    if (k != c_inst_args.kmer_length) {
        c_err << "ERROR: KMC database and provided k-mer lengths are different." << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    for (std::size_t it = 0; it < c_inst_args.read_files_data.size(); ++it) {
        C_correct_errors reads_corrector(c_inst_args, c_inst_args.read_files_data[it], check_inputs[it], c_inst_args.accept_filtered_with_long_kmers);

        c_log << c_inst_args.read_files_data[it].input_name << std::endl;

        Timer correct_errors_in_reads_timer;
        correct_errors_in_reads_timer.startTimer();

        // correct errors in reads
        if (c_inst_args.use_long_kmer) {
            reads_corrector.correct_errors_in_a_file(kmc_file, &kmc_long_file);
        }
        else {
            reads_corrector.correct_errors_in_a_file(kmc_file, nullptr);
        }

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

    if (!c_inst_args.reuse_kmc_db) {
        PrepareKMCDb prepareKMCDb(c_inst_args);
        prepareKMCDb.removeFilteredKMCDatabase();
    }
}