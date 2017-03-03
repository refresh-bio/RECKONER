/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 1.0
*
*/

#include "DetermineParameters.hpp"
#include "FileReader.h"
#include "Log.h"
#include "HistogramAnalyzer.h"
#include <string>
#include <cmath>
#include <algorithm>



inline double DetermineParameters::phredToProb(char symbol, std::size_t _qualityScoreOffset) {
    symbol -= static_cast<char>(_qualityScoreOffset);
    return std::pow(static_cast<double>(10.0), static_cast<double>(-symbol / 10.0));
}

void DetermineParameters::findReadsStats(double& errorRate, double& readLength) {
    double qualitySum = 0.0;
    unsigned qualitiesNumber = 0;
    unsigned readsNumber = 0;

    for (unsigned i = 0; i < c_inst_args.read_files_names.size(); ++i) {
        FileReader currentFile;
        currentFile.setFileName(c_inst_args.read_files_names[i], c_inst_args.read_files_types[i]);
        if (!currentFile.openFile(FileReader::READ)) {
            std::cerr << "ERROR: Cannot open " << c_inst_args.read_files_names[i] << " to determine the k-mer length." << std::endl;
            Log::get_stream() << "ERROR: Cannot open " << c_inst_args.read_files_names[i] << " to determine the k-mer length." << std::endl;
            exit(EXIT_FAILURE);
        }

        std::string line;
        unsigned readsTested = 0;

        // header
        currentFile.getLine(line);

        while (!currentFile.eof()) {
            // sequence
            currentFile.getLine(line);
            // +
            currentFile.getLine(line);
            // quality score
            currentFile.getLine(line);

            for (unsigned j = 0; j < line.length(); ++j) {
                qualitySum += phredToProb(line[j], check_inputs[i].quality_score_offset);
            }
            qualitiesNumber += static_cast<unsigned>(line.length());

            ++readsTested;
            if (readsTested >= DETERMINE_READS_FOR_QUALITY) {
                break;
            }

            // header
            currentFile.getLine(line);
        }

        readsNumber += readsTested;
    }

    errorRate = qualitySum / qualitiesNumber;
    readLength = static_cast<double>(qualitiesNumber) / readsNumber;
}

std::size_t DetermineParameters::determineKmerLength() {
    double errorRate = 0;
    double readLength = 0;

    findReadsStats(errorRate, readLength);

    errorRate *= 100.0;

    double a = 0.8;
    double b = 9 + errorRate;
    double c = 2 + static_cast<double>(readLength) / 100;

    double kMin = 20;
    double kMax = 0.2 * readLength + 30;

    double k = (a * log2(c_inst_args.genome_size) - b) * c;

    k = std::min(k, kMax);
    k = std::max(k, kMin);

    return static_cast<std::size_t>(std::round(k));
}

std::size_t DetermineParameters::determineGenomeSize() {
    char* kmcBuffer = new char[KMC_STDOUT_BUFFER_SIZE];
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory);
    if (!runExternal.runKMCToBuffer(PROBE_KMER_LENGTH, c_inst_args.read_files_names, c_inst_args.kmc_determine_params_database_name, c_inst_args.prefix, kmcBuffer, KMC_STDOUT_BUFFER_SIZE)) {
        std::cerr << "ERROR: failed to count k-mers." << std::endl;
        Log::get_stream() << "ERROR: failed to count k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::size_t allKmers;
    std::size_t allReads;

    if (!parseKMCOutput(allKmers, allReads, kmcBuffer, KMC_STDOUT_BUFFER_SIZE)) {
        std::cerr << "ERROR: cannot obtain k-mers statistics." << std::endl;
        Log::get_stream() << "ERROR: cannot obtain k-mers statistics." << std::endl;
        exit(EXIT_FAILURE);
    }
    delete[] kmcBuffer;

    unsigned histoMaximum;
    HistogramAnalyzer histogramAnalyzer;
    histoMaximum = histogramAnalyzer.getHistogramPeak(c_inst_args.kmc_determine_params_database_name);

    removeKMCDatabase();

    unsigned readLength = static_cast<unsigned>(std::round((allKmers / static_cast<double>(allReads + PROBE_KMER_LENGTH - 1)) + PROBE_KMER_LENGTH - 1));
    double coverage = static_cast<float>(histoMaximum * readLength) / static_cast<double>(readLength - PROBE_KMER_LENGTH + 1);

    //std::cout << "All kmers   " << allKmers << std::endl;
    //std::cout << "All reads   " << allReads << std::endl;
    //std::cout << "Maximum     " << histoMaximum << std::endl;
    //std::cout << "Read length " << readLength << std::endl;
    //std::cout << "Coverage    " << coverage << std::endl;

    std::size_t genomeSize = static_cast<std::size_t>(allReads * readLength / coverage);

    return genomeSize;
}

bool DetermineParameters::getValueFromString(const std::string& line, const std::string& prefixText, std::size_t& resultValue) {
    std::size_t pos = line.find(prefixText);
    if (pos == std::string::npos) {
        return false;
    }
    pos = line.rfind(":");

    std::istringstream sstream(line.substr(pos + 1));
    std::size_t result;
    sstream >> result;
    if (sstream) {
        resultValue = result;
        return true;
    }
    return false;
}

bool DetermineParameters::parseKMCOutput(std::size_t& resultAllKmers, std::size_t& resultAllReads, char* buffer, const int bufferSize) {
    bool allKmersSet = false;
    bool allReadsSet = false;
    std::size_t allReads, allKmers;

    std::istringstream sstream(buffer);
    std::string line;
    while (getline(sstream, line)) {
        if (getValueFromString(line, "Total no. of k-mers", allKmers)) {
            allKmersSet = true;
        }
        else if (getValueFromString(line, "Total no. of reads", allReads)) {
            allReadsSet = true;
        }
    }


    if (allKmersSet && allReadsSet) {
        resultAllKmers = allKmers;
        resultAllReads = allReads;
        return true;
    }

    return false;
}

void DetermineParameters::removeKMCDatabase() {
    RunExternal::removeKMCFiles(c_inst_args.kmc_determine_params_database_name);
}
