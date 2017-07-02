/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 1.1
*
*/

#include "RunExternal.h"
#include "PrepareKMCDb.h"
#include "Log.h"
#include "HistogramAnalyzer.h"
#include <iostream>



void PrepareKMCDb::countKmers() {
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory);

    if (!runExternal.runKMC(static_cast<int>(c_inst_args.kmer_length), c_inst_args.read_files_names, c_inst_args.kmc_database_name, c_inst_args.kmc_list_file_name, c_inst_args.prefix)) {
        std::cerr << "ERROR: failed to count k-mers." << std::endl;
        Log::get_stream() << "ERROR: failed to count k-mers." << std::endl;
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
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory);

    if (!runExternal.runKMCTools(cutoff, c_inst_args.kmc_database_name, c_inst_args.kmc_filtered_database_name)) {
        std::cerr << "ERROR: failed to remove untrusted k-mers." << std::endl;
        Log::get_stream() << "ERROR: failed to remove untrusted k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void PrepareKMCDb::removeKMCDatabase() {
    RunExternal::removeKMCFiles(c_inst_args.kmc_database_name);
}

void PrepareKMCDb::removeFilteredKMCDatabase() {
    RunExternal::removeKMCFiles(c_inst_args.kmc_filtered_database_name);
}
