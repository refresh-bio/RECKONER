/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 2.1
*
*/

#include "HistogramAnalyzer.h"
#include "define.hpp"



bool HistogramAnalyzer::buildHistogram(const std::string& kmcDatabase, const unsigned histogramSize, std::vector<uint64>& histogram) {
    CKMCFile kmcFile;
    if (!kmcFile.OpenForListing(kmcDatabase)) {
        return false;
    }

    CKMCFileInfo kmcFileInfo;
    if (!kmcFile.Info(kmcFileInfo)) {
        return false;
    }

    CKmerAPI kmer(kmcFileInfo.kmer_length);
    uint32 count;
    unsigned iCount;

    histogram.clear();
    histogram.resize(histogramSize, 0);

    while (kmcFile.ReadNextKmer(kmer, count)) {
        iCount = static_cast<int> (std::round(count));
        if (iCount >= histogramSize) {
            iCount = histogramSize - 1;
        }

        ++histogram[iCount];
    }

    return true;
}

unsigned HistogramAnalyzer::getCutoffThreshold(const std::string& kmcDatabase) {
    std::vector<uint64> histogram;

    if (!buildHistogram(kmcDatabase, HISTOGRAM_SIZE, histogram)) {
        c_err << "ERROR: cannot open KMC files to determine cutoff threshold." << std::endl;
        exit(EXIT_FAILURE);
    }

    uint64 firstMin = histogram[2];

    unsigned cutoff = 2;
    for (unsigned i = 3; i < HISTOGRAM_SIZE; ++i) {
        if (firstMin >= histogram[i]) {
            firstMin = histogram[i];
        }
        else {
            cutoff = i - 1;
            break;
        }
    }

    return (cutoff > MAX_CUTOFF ? MAX_CUTOFF : cutoff);
}

unsigned HistogramAnalyzer::getHistogramPeak(const std::string& kmcDatabase) {
    std::vector<uint64> histogram;

    if (!buildHistogram(kmcDatabase, HISTOGRAM_SIZE, histogram)) {
        c_err << "ERROR: cannot open KMC files to determine k-mers histogram peak." << std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned firstMin = 2;
    for (unsigned i = 3; i < HISTOGRAM_SIZE; ++i) {
        if (histogram[firstMin] >= histogram[i]) {
            firstMin = i;
        }
        else {
            break;
        }
    }

    unsigned maximum = firstMin + 1;
    for (unsigned i = maximum; i < HISTOGRAM_SIZE; ++i) {
        if (histogram[maximum] < histogram[i]) {
            maximum = i;
        }
    }

    return maximum;
}

