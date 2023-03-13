/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 2.0
*
*/

#include "RunExternal.h"
#include "PrepareKMCDb.h"
#include "Log.h"
#include "HistogramAnalyzer.h"
#include <iostream>



void PrepareKMCDb::countKmers() {
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory, c_inst_args.kmc_ram);

    // extract input files names
    std::vector<std::string> read_files_names;
    for (auto& it : c_inst_args.read_files_data) {
        read_files_names.push_back(it.input_name);
    }

    if (!runExternal.runKMC(static_cast<int>(c_inst_args.kmer_length), read_files_names, c_inst_args.kmc_database_name, c_inst_args.kmc_list_file_name, c_inst_args.prefix)) {
        c_err << "ERROR: failed to count k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void PrepareKMCDb::countLongKmers() {
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory, c_inst_args.kmc_ram);

    // extract input files names
    std::vector<std::string> read_files_names;
    for (auto& it : c_inst_args.read_files_data) {
        read_files_names.push_back(it.input_name);
    }

    if (!runExternal.runKMC(static_cast<int>(c_inst_args.long_kmer_length), read_files_names, c_inst_args.kmc_long_database_name, c_inst_args.kmc_list_file_name, c_inst_args.prefix, 1)) {
        c_err << "ERROR: failed to count long k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }
}

unsigned PrepareKMCDb::determineCutoff() {
    unsigned cutoff;

    HistogramAnalyzer histogramAnalyzer;
    cutoff = histogramAnalyzer.getCutoffThreshold(c_inst_args.kmc_database_name);

    return cutoff;
}

void PrepareKMCDb::filterKmers(unsigned cutoff) {
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory, c_inst_args.kmc_ram);

    if (!runExternal.runKMCTools(cutoff, c_inst_args.kmc_database_name, c_inst_args.kmc_filtered_database_name)) {
        c_err << "ERROR: failed to remove untrusted k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void PrepareKMCDb::removeKMCDatabase() {
    RunExternal::removeKMCFiles(c_inst_args.kmc_database_name);
}

void PrepareKMCDb::removeFilteredKMCDatabase() {
    RunExternal::removeKMCFiles(c_inst_args.kmc_filtered_database_name);
    if (c_inst_args.use_long_kmer) {
        RunExternal::removeKMCFiles(c_inst_args.kmc_long_database_name);
    }
}
