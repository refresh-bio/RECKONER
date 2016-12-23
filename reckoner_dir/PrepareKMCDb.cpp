/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* 0.2.1
*
*/

#include "RunExternal.h"
#include "PrepareKMCDb.h"
#include "Log.h"
#include "HistogramAnalyzer.h"
#include <iostream>

void PrepareKMCDb::run(const C_arg& c_inst_args) {
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory);

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "K - MER COUNTING" << std::endl;
    std::cout << "##################################################" << std::endl;
    Log::get_stream() << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    Log::get_stream() << "K - MER COUNTING" << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    if (!runExternal.runKMC(c_inst_args.kmer_length, c_inst_args.read_files_names, c_inst_args.read_files_names.front() + ".lst.tmp", c_inst_args.prefix)) {
        std::cerr << "ERROR: failed to count k-mers." << std::endl;
        Log::get_stream() << "ERROR: failed to count k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "DETERMINING CUTOFF THRESHOLD" << std::endl;
    std::cout << "##################################################" << std::endl;
    Log::get_stream() << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    Log::get_stream() << "DETERMINING CUTOFF THRESHOLD" << std::endl;
    Log::get_stream() << "##################################################" << std::endl;

    unsigned cutoff;
    HistogramAnalyzer histogramAnalyzer;
    cutoff = histogramAnalyzer.getCutoffThreshold(c_inst_args.read_files_names.front() + ".lst.tmp");

    std::cout << "Cutoff: " << cutoff << std::endl;

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "REMOVING UNTRUSTED K - MERS" << std::endl;
    std::cout << "##################################################" << std::endl;
    Log::get_stream() << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    Log::get_stream() << "REMOVING UNTRUSTED K - MERS" << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    if (!runExternal.runKMCTools(cutoff, c_inst_args.read_files_names.front() + ".lst.tmp", c_inst_args.read_files_names.front() + ".lst")) {
        std::cerr << "ERROR: failed to remove untrusted k-mers." << std::endl;
        Log::get_stream() << "ERROR: failed to remove untrusted k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << std::endl;
}
