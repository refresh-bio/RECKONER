/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.0
 * 
 */

#ifndef _DEFINE_H
#define _DEFINE_H

// C++ libraries
#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <iterator>
#include <cctype>


#include <sys/stat.h>
#include <sys/types.h>



// definitions
#define VERSION                    "1.0"
#define NUM_NEOCLEOTIDE            4
#define A                          0
#define C                          1
#define G                          2
#define T                          3
#define MIN_KMER_LENGTH            3
#define PHRED33                    33
#define PHRED64                    64
#define MIN_33_SCORE               33
#define MAX_33_SCORE               74
#define MIN_64_SCORE               59 // for qualities <-5, 40>
#define MAX_64_SCORE               104
#define MAX_DETECT_SCORE_DIFF      5
#define MIN_SCORE                  33
#define MIN_NON_SOLID_LENGTH       2
#define MIN_SOLID_LENGTH           2
#define MIN_QS_DIFF                10
#define MAX_EXTENSION              5
#define QS_CUTOFF                  10
#define MAX_LOW_QS_BASES           4
#define FP_SUSPECT_LENGTH          1
#define SOLID_REGION_ADJUST_RANGE  4
#define SUBST_CHAR                 'A'
#define MAX_ERROR_RATE             0.5
#define MIN_BEST_KMER_QUALITY      0.0
#define COVERING_KMERS_WEIGHT      1.0
#define EXTENSION_KMERS_WEIGHT     0.5
#define MAX_CHANGES_IN_REGION_RATIO 0.5
//#define LIMIT_MODIFICATIONS
#define MAX_LOW_QS_INDEXES_COMB    6
#define MAX_FIRST_KMER_POSSIBILITIES 5
#define READ_LINES                 4
#define MAX_EXTEND_CORRECTION_PATHS 100
#define MAX_FIRST_KMER_CORRECTION_PATHS 30
#define MAX_N_RATIO                 0.3
#define CHECK_MAX_CHANGES          10000
#define MAX_CHECK_FIRST_KMER_NESTING (MAX_LOW_QS_BASES > MAX_LOW_QS_INDEXES_COMB ? MAX_LOW_QS_BASES : MAX_LOW_QS_INDEXES_COMB)
#define LIST_FILE_EXTENSION        ".lst"
#define TEMP_EXTENSION             ".tmp"
#define DETERMINE_PARAMS_EXTENSION ".gen"
#define LOG_FILE_EXTENSION         ".log"
#define CORRECTION_FILE_EXTENSION  ".error-correction"
#define CORRECTED_FILE_EXTENSION   ".corrected"
#define GZ_FILE_EXTENSION          ".gz"
#define DEFAULT_MAX_KMC_MEMORY_USAGE 4
#define DETERMINE_READS_FOR_QUALITY 100000
#define HISTOGRAM_SIZE              256
#define MAX_CUTOFF                  5
#define PART_SIZE                   (8 << 20)
#define PART_BUFFERS_PER_THREAD     2
#ifdef WIN32 // Windows
#define KMC_EXECUTABLE_NAME "kmc.exe"
#define KMC_TOOLS_EXECUTABLE_NAME "kmc_tools.exe"
#else // Linux
#define KMC_EXECUTABLE_NAME "kmc"
#define KMC_TOOLS_EXECUTABLE_NAME "kmc_tools"
#endif

static const char NEOCLEOTIDE[NUM_NEOCLEOTIDE] = {'A', 'C', 'G', 'T'};

#endif
