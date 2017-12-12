/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 1.1.1
*
*/

#ifndef DETERMINE_PARAMETERS_HPP
#define DETERMINE_PARAMETERS_HPP

#include "RunExternal.h"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include <vector>
#include <sstream>



class DetermineParameters {
private:
    static const int PROBE_KMER_LENGTH = 25;
    static const int KMC_STDOUT_BUFFER_SIZE = 50000;

    double phredToProb(char symbol, std::size_t _qualityScoreOffset);
    void findReadsStats(double& errorRate, double& readLength);
    bool getValueFromString(const std::string& line, const std::string& prefixText, std::size_t& resultValue);
    bool parseKMCOutput(std::size_t& resultAllKmers, std::size_t& resultAllReads, char* buffer, const int bufferSize);
    void removeKMCDatabase();

    const C_arg& c_inst_args;
    const std::vector<C_check_read>& check_inputs;

public:
    DetermineParameters(const C_arg& _c_inst_args, const std::vector<C_check_read>& _check_inputs) :
        c_inst_args(_c_inst_args), check_inputs(_check_inputs) {}

    std::size_t determineKmerLength();
    std::size_t determineGenomeSize();
};


#endif
