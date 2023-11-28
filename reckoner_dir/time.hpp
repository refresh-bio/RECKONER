/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.1
 * 
 */

#ifndef _TIME_LOCAL_H
#define _TIME_LOCAL_H

#define _CRT_SECURE_NO_WARNINGS

#include "define.hpp"
#include <ctime>
#include <vector>
#include <string>



#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4996)
#endif
class Timer {
private:
    std::time_t startTime;
    std::time_t stopTime;

public:
    Timer() : startTime(0), stopTime(0) {}

    void startTimer() {
        time(&startTime);
    }
    void stopTimer() {
        time(&stopTime);
    }

    std::string getStartTime() {
        std::string result = ctime(&startTime);
        result.pop_back();
        return result;
    }
    std::string getStopTime() {
        std::string result = ctime(&stopTime);
        result.pop_back();
        return result;
    }

    bool wasTimeMeasured() {
        return startTime != 0 && stopTime != 0;
    }
};
#ifdef _MSC_VER
#pragma warning(pop)
#endif




//----------------------------------------------------------------------
// C_time
//----------------------------------------------------------------------

class C_time {
public:
    // variables

    // check read files
    std::vector<Timer> vector_check_read_file;

    // determine genome size and k-mer length
    Timer determine_parameters;

    // count k-mers
    Timer kmer_count;

    // count long k-mers
    Timer long_kmer_count;

    // determine cutoff threshold
    Timer determine_cutoff_threshold;

    // remove untrusted k-mers
    Timer remove_untrusted_kmers;

    // correct errors in reads
    std::vector<Timer> vector_correct_errors_in_reads;
};



#endif
