/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* 0.2.1
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
    symbol -= _qualityScoreOffset;
    return std::pow(static_cast<double>(10.0), static_cast<double>(-symbol / 10.0));
}

double DetermineParameters::finedMeanErrorRate() {
    double qualitySum = 0.0;
    unsigned qualitiesNumber = 0;

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
            qualitiesNumber += line.length();

            ++readsTested;
            if (readsTested >= DETERMINE_READS_FOR_QUALITY) {
                break;
            }

            // header
            currentFile.getLine(line);
        }
    }

    return qualitySum / qualitiesNumber;
}

std::size_t DetermineParameters::determineKmerLength() {
    double errorRate = finedMeanErrorRate();
    errorRate *= 100.0;

    double minK = 29.0 - 2.0 * errorRate;

    double a = std::max(6.0 - errorRate, 1.0);
    double b = 55.0 - 10.0 * errorRate;

    double k = std::max(a * std::log2(c_inst_args.genome_size / 1000.0) - b, minK);

    return static_cast<std::size_t>(std::round(k));
}

std::size_t DetermineParameters::determineGenomeSize() {
    char* kmcBuffer = new char[KMC_STDOUT_BUFFER_SIZE];
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory);
    if (!runExternal.runKMCToBuffer(PROBE_KMER_LENGTH, c_inst_args.read_files_names, c_inst_args.read_files_names.front() + ".gen", c_inst_args.prefix, kmcBuffer, KMC_STDOUT_BUFFER_SIZE)) {
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
    histoMaximum = histogramAnalyzer.getHistogramPeak(c_inst_args.read_files_names.front() + ".gen");

    unsigned readLength = std::round((allKmers / static_cast<double>(allReads + PROBE_KMER_LENGTH - 1)) + PROBE_KMER_LENGTH - 1);
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

void DetermineParameters::perform_determine_parameters() {
    if (!c_inst_args.is_kmer_length_user_defined) {
        std::cout << std::endl;
        std::cout << "##################################################" << std::endl;
        std::cout << "DETERMINING PARAMETERS" << std::endl;
        std::cout << "##################################################" << std::endl;
        Log::get_stream() << std::endl;
        Log::get_stream() << "##################################################" << std::endl;
        Log::get_stream() << "DETERMINING PARAMETERS" << std::endl;
        Log::get_stream() << "##################################################" << std::endl;

        if (!c_inst_args.is_genome_size_user_defined) {
            c_inst_args.genome_size = determineGenomeSize();
            std::cout << "Estimated genome size: " << c_inst_args.genome_size << std::endl;
        }
        c_inst_args.kmer_length = determineKmerLength();
        std::cout << "Determined k-mer length: " << c_inst_args.kmer_length << std::endl;
    }
}
