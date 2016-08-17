/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 0.2
*
*/

#include "RunExternal.h"
#include <cstring>
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <iterator>
#include <algorithm>
#include <sstream>
#endif


bool RunExternal::runKMCTools(int cutoff, const std::string& inputFileName, const std::string& outputFileName) {
    std::string args;

    args += "reduce";
    args += " " + inputFileName;

    std::ostringstream sstream;
    sstream << cutoff;
    args += " -ci" + sstream.str();

    args += " " + outputFileName;

    unsigned long result;
    if (!runCommand(KMC_TOOLS_EXECUTABLE_NAME, args, result) || result != 0) {
        std::cerr << "ERROR: cannot run kmc_tools." << std::endl;
        Log::get_stream() << "ERROR: cannot run kmc_tools." << std::endl;
        return false;
    }

    return true;
}

bool RunExternal::runKMC(int kmerLength, const std::vector<std::string>& inputFilesNames, const std::string& outputFileName, const std::string& tempName) {
    if (!createDirectory(tempName)) {
        std::cerr << "ERROR: cannot create KMC temp directory." << std::endl;
        Log::get_stream() << "ERROR: cannot create KMC temp directory." << std::endl;
        return false;
    }

    std::string listFileName = inputFilesNames.front();
    listFileName += ".lst";

    std::ofstream listFileStream(listFileName);
    if (!listFileStream.is_open()) {
        std::cerr << "ERROR: cannot create list of KMC input files." << std::endl;
        Log::get_stream() << "ERROR: cannot create list of KMC input files." << std::endl;
        return false;
    }

    for (const std::string& inputFileName : inputFilesNames) {
        listFileStream << inputFileName << "\n";
    }
    listFileStream.close();

    std::string args;

    std::ostringstream sstream;
    sstream << kmerLength;
    args += "-k" + sstream.str();

    sstream.str("");
    sstream << maxMemory;
    args += " -m" + sstream.str();

    sstream.str("");
    sstream << nThreads;
    args += " -t" + sstream.str();

    args += " -ci2";
    args += " @" + listFileName;
    args += " " + outputFileName;
    args += " " + tempName;

    unsigned long result;
    if (!runCommand(KMC_EXECUTABLE_NAME, args, result) || result != 0) {
        std::cerr << "ERROR: cannot run KMC." << std::endl;
        Log::get_stream() << "ERROR: cannot run KMC." << std::endl;
        return false;
    }

    std::remove(listFileName.c_str());

    return true;
}

bool RunExternal::runKMCToBuffer(std::size_t kmerLength, const std::vector<std::string>& inputFilesNames, const std::string& outputFileName, const std::string& tempName, char* buffer, const int bufferSize) {
    if (!createDirectory(tempName)) {
        std::cerr << "ERROR: cannot create KMC temp directory." << std::endl;
        Log::get_stream() << "ERROR: cannot create KMC temp directory." << std::endl;
        return false;
    }

    std::string listFileName = inputFilesNames.front();
    listFileName += ".lst";

    std::ofstream listFileStream(listFileName);
    if (!listFileStream.is_open()) {
        std::cerr << "ERROR: cannot create list of KMC input files." << std::endl;
        Log::get_stream() << "ERROR: cannot create list of KMC input files." << std::endl;
        return false;
    }

    for (const std::string& inputFileName : inputFilesNames) {
        listFileStream << inputFileName << "\n";
    }
    listFileStream.close();

    std::string args;

    std::ostringstream sstream;
    sstream << kmerLength;
    args += "-k" + sstream.str();

    sstream.str("");
    sstream << maxMemory;
    args += " -m" + sstream.str();

    sstream.str("");
    sstream << nThreads;
    args += " -t" + sstream.str();

    args += " -ci2";
    args += " @" + listFileName;
    args += " " + outputFileName;
    args += " " + tempName;

    unsigned long result;
    if (!runCommand(KMC_EXECUTABLE_NAME, args, result, buffer, bufferSize) || result != 0) {
        std::cerr << "ERROR: cannot run KMC." << std::endl;
        Log::get_stream() << "ERROR: cannot run KMC." << std::endl;
        return false;
    }

    std::remove(listFileName.c_str());

    return true;
}

void RunExternal::removeKMCFiles(const std::string& fileName) {
    std::remove((fileName + ".kmc_pre").c_str());
    std::remove((fileName + ".kmc_suf").c_str());
}

#ifdef WIN32 //Windows

bool RunExternal::runCommand(const std::string& command, const std::string& args, unsigned long& processResult, char* buffer /*= NULL*/, const int bufferSize /*= 0*/) {
    PROCESS_INFORMATION processInfo;
    STARTUPINFO processStartupInfo;

    ZeroMemory(&processInfo, sizeof(PROCESS_INFORMATION));

    ZeroMemory(&processStartupInfo, sizeof(STARTUPINFO));
    processStartupInfo.cb = sizeof(STARTUPINFO);
    processStartupInfo.hStdError = GetStdHandle(STD_OUTPUT_HANDLE);

    HANDLE pipeReadHandle = NULL;
    HANDLE pipeWriteHandle = NULL;

    if (buffer != NULL){
        if (!createPipe(pipeReadHandle, pipeWriteHandle)) {
            return false;
        }
        processStartupInfo.hStdOutput = pipeWriteHandle;
    }
    else {
        processStartupInfo.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    }

    processStartupInfo.hStdInput = NULL;
    processStartupInfo.dwFlags |= STARTF_USESTDHANDLES;

    // Fullfill the ANSI specifications by passing command name as argv[0].
    char* cArgs = new char[command.length() + 1 + args.length() + 1];
    std::size_t copied = command.copy(cArgs, command.length());
    cArgs[copied] = ' ';
    ++copied;
    copied += args.copy(cArgs + copied, args.length());
    cArgs[copied] = '\0';

    BOOL createProcessSuccess = CreateProcess(command.c_str(), cArgs, NULL, NULL, TRUE, 0, NULL, NULL, &processStartupInfo, &processInfo);

    delete[] cArgs;

    if (!createProcessSuccess) {
        return false;
    }
    else {
        WaitForSingleObject(processInfo.hProcess, INFINITE);

        if (buffer != NULL) {
            DWORD bytesRead;

            BOOL pipeReadSuccess = PeekNamedPipe(pipeReadHandle, buffer, bufferSize - 1, &bytesRead, NULL, NULL);

            if (!pipeReadSuccess || bytesRead == 0) {
                return false;
            }
            buffer[bytesRead] = '\0';

            destroyPipe(pipeReadHandle, pipeWriteHandle);
        }

        if (!GetExitCodeProcess(processInfo.hProcess, &processResult)) {
            return false;
        }

        CloseHandle(processInfo.hProcess);
        CloseHandle(processInfo.hThread);

        return true;
    }
}

bool RunExternal::createPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle) {
    SECURITY_ATTRIBUTES securityAttributes;
    securityAttributes.nLength = sizeof(SECURITY_ATTRIBUTES);
    securityAttributes.bInheritHandle = TRUE;
    securityAttributes.lpSecurityDescriptor = NULL;

    if (!CreatePipe(&pipeReadHandle, &pipeWriteHandle, &securityAttributes, 0)) {
        return false;
    }

    if (!SetHandleInformation(pipeReadHandle, HANDLE_FLAG_INHERIT, 0)) {
        return false;
    }

    return true;
}

void RunExternal::destroyPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle) {
    CloseHandle(pipeReadHandle);
    CloseHandle(pipeWriteHandle);

    pipeReadHandle = NULL;
    pipeWriteHandle = NULL;
}

bool RunExternal::createDirectory(const std::string& dirName) {
    if (_mkdir(dirName.c_str()) != 0) {
        if (errno != EEXIST) {
            return false;
        }
    }

    return true;
}

#else // Linux

bool RunExternal::runCommand(const std::string& command, const std::string& args, unsigned long& processResult, char* buffer /*= NULL*/, const int bufferSize /*= 0*/) {
    if (buffer == NULL) {
        int pid = fork();
        if (pid == -1) {
            return false;
        }
        else if (pid == 0) {
            // Fullfill the ANSI specifications by passing command name as argv[0].
            std::istringstream sstream(command + " " + args);
            std::vector<std::string> splittedArgs;
            copy(std::istream_iterator<std::string>(sstream), std::istream_iterator<std::string>(), back_inserter(splittedArgs));
            char** cSplittedArgs = new char*[splittedArgs.size() + 1];
            for (unsigned i = 0; i < splittedArgs.size(); ++i) {
                cSplittedArgs[i] = new char[splittedArgs[i].length() + 1];
                std::size_t copied = splittedArgs[i].copy(cSplittedArgs[i], splittedArgs[i].length());
                cSplittedArgs[i][copied] = '\0';
            }
            cSplittedArgs[splittedArgs.size()] = NULL;

            execv(cSplittedArgs[0], cSplittedArgs);
            for (unsigned i = 0; i < splittedArgs.size(); ++i) {
                delete[] cSplittedArgs[i];
            }
            delete[] cSplittedArgs;
            exit(EXIT_FAILURE);
        }
        else {
            int status;
            if (waitpid(pid, &status, 0) < 0) {
                return false;
            }

            processResult = WEXITSTATUS(status);
            return true;
        }
    }
    else {
        std::string childCommand = command + " " + args;

        FILE* pipe = popen(childCommand.c_str(), "r");
        if (pipe == NULL) {
            return false;
        }

        std::size_t bytesRead = fread(buffer, sizeof(char), bufferSize - 1, pipe);
        if (bytesRead == 0) {
            return false;
        }
        buffer[bytesRead] = '\0';

        int result = pclose(pipe);
        processResult = WEXITSTATUS(result);

        return true;
    }
}

bool RunExternal::createDirectory(const std::string& dirName) {
    if (mkdir(dirName.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) != 0) {
        if (errno != EEXIST) {
            return false;
        }
    }

    return true;
}

#endif
