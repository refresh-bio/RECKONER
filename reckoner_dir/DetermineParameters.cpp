/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 2.1
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

    for (unsigned i = 0; i < c_inst_args.read_files_data.size(); ++i) {
        FileReader currentFile;
        currentFile.setFileName(c_inst_args.read_files_data[i].input_name, c_inst_args.read_files_data[i].type);
        if (!currentFile.openFile(FileReader::READ)) {
            c_err << "ERROR: Cannot open " << c_inst_args.read_files_data[i].input_name << " to determine the k-mer length." << std::endl;
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
    double k = K_LEN_A * log2(c_inst_args.genome_size) + K_LEN_B;

    k = std::max(k, static_cast<double>(K_MIN));

    return static_cast<std::size_t>(std::round(k));
}

std::size_t DetermineParameters::determineGenomeSize() {
    char* kmcBuffer = new char[KMC_STDOUT_BUFFER_SIZE];
    RunExternal runExternal(c_inst_args.n_threads, c_inst_args.kmc_memory, c_inst_args.kmc_ram);

    // extract input files names
    std::vector<std::string> read_files_names;
    for (auto& it : c_inst_args.read_files_data) {
        read_files_names.push_back(it.input_name);
    }

    if (!runExternal.runKMCToBuffer(PROBE_KMER_LENGTH, read_files_names, c_inst_args.kmc_determine_params_database_name, c_inst_args.kmc_list_file_name, c_inst_args.prefix, kmcBuffer, KMC_STDOUT_BUFFER_SIZE)) {
        c_err << "ERROR: failed to count k-mers." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::size_t allKmers;
    std::size_t allReads;

    if (!parseKMCOutput(allKmers, allReads, kmcBuffer, KMC_STDOUT_BUFFER_SIZE)) {
        c_err << "ERROR: cannot obtain k-mers statistics." << std::endl;
        exit(EXIT_FAILURE);
    }
    delete[] kmcBuffer;

    unsigned histoMaximum;
    HistogramAnalyzer histogramAnalyzer;
    histoMaximum = histogramAnalyzer.getHistogramPeak(c_inst_args.kmc_determine_params_database_name);

    removeKMCDatabase();

    unsigned readLength = static_cast<unsigned>(std::round((allKmers / static_cast<double>(allReads + PROBE_KMER_LENGTH - 1)) + PROBE_KMER_LENGTH - 1));
    double coverage = static_cast<double>(histoMaximum * readLength) / static_cast<double>(readLength - PROBE_KMER_LENGTH + 1);

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
    std::size_t allReads = 0, allKmers = 0;

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
