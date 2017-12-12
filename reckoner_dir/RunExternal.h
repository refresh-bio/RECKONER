/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 1.1.1
*
*/

#ifndef RUNEXTERNAL_H
#define	RUNEXTERNAL_H

#ifdef WIN32
#define NOMINMAX // to use std::max
#include <Windows.h>
#else
#include <sys/stat.h>
#endif

#include "Log.h"
#include "define.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>



class RunExternal {
    std::size_t nThreads;
    std::size_t maxMemory;
    static const int STDOUT_BUFFER_SIZE = 100;

public:
    RunExternal(std::size_t _nThreads, std::size_t _maxMemory) : nThreads(_nThreads), maxMemory(_maxMemory) {}

    bool runKMCTools(int cutoff, const std::string& inputFileName, const std::string& outputFileName);
    bool runKMC(int kmerLength, const std::vector<std::string>& inputFilesNames, const std::string& outputFileName, const std::string& listFileName, const std::string& tempName);
    bool runKMCToBuffer(std::size_t kmerLength, const std::vector<std::string>& inputFilesNames, const std::string& outputFileName, const std::string& listFileName, const std::string& tempName, char* buffer, const int bufferSize);

    static void removeKMCFiles(const std::string& fileName);

private:
#ifdef WIN32 // Windows

    // buffer - output for child process' stdout and stderr, if NULL, the result is sent to stdout
    bool runCommand(std::string command, const std::string& args, unsigned long& processResult, char* buffer = NULL, const int bufferSize = 0);

    bool createPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle);
    void destroyPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle);

public:
    static bool createDirectory(const std::string& dirName);

#else // Linux

    // buffer - output for child process' stdout and stderr, if NULL, the result is sent to stdout
    bool runCommand(std::string command, const std::string& args, unsigned long& processResult, char* buffer = NULL, const int bufferSize = 0);

public:
    static bool createDirectory(const std::string& dirName);

#endif
};

#endif	/* RUNEXTERNAL_H */
