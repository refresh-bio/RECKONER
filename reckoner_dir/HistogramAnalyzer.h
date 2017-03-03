/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 1.0
*
*/

#ifndef HISTOGRAMANALYZER_H
#define HISTOGRAMANALYZER_H

#include <kmc_api/kmer_api.h>
#include <kmc_api/kmc_file.h>
#include <string>
#include <vector>
#include <cmath>



class HistogramAnalyzer {
private:
    bool buildHistogram(const std::string& kmcDatabase, const unsigned histogramSize, std::vector<uint64>& histogram);

public:
    unsigned getCutoffThreshold(const std::string& kmcDatabase);
    unsigned getHistogramPeak(const std::string& kmcDatabase);
};

#endif /* HISTOGRAMANALYZER_H */
