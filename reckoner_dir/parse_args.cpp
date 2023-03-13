/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 2.0
 * 
 */

#include "RunExternal.h"
#undef T // workaround for newer g++s issue
#include <thread>
#include "parse_args.hpp"
#include "Log.h"
#include <string.h>
#include <algorithm>



//----------------------------------------------------------------------
// Constructor.
//----------------------------------------------------------------------

C_arg::C_arg(int argc, char** argv) :
    is_kmer_length_user_defined(false),
    kmer_length(0),
    long_kmer_length(0),
    is_genome_size_user_defined(false),
    genome_size(0),
    extend(MAX_EXTENSION),
    n_threads(std::thread::hardware_concurrency()),
    kmc_memory(DEFAULT_MAX_KMC_MEMORY_USAGE),
    kmc_ram(false),
    reuse_kmc_db(false),
    use_long_kmer(false),
    long_kmer_ratio(DEFAULT_LONG_KMER_RATIO),
    accept_filtered_with_long_kmers(false),
    correct_indels(true),
    change_headers_length_tag(false),
    mark_corrected(false),
    verbose(false),
    num_args(argc),
    args(argv) {
    if (n_threads == 0) {
        n_threads = 1;
    }
    read_args();
}



//----------------------------------------------------------------------
// Extracts extension from file name and detects file type.
//----------------------------------------------------------------------

ReadFileData::ReadFileData(const std::string& _input_name)
    : input_name(_input_name) {
    std::size_t slash_pos = _input_name.find_last_of("/\\");
    std::string result_file_name(slash_pos == std::string::npos ? _input_name : _input_name.substr(slash_pos + 1));
    std::string last_extension, pre_last_extension;

    std::size_t dot_pos = result_file_name.find_last_of(".");
    if (dot_pos != std::string::npos) {
        last_extension = result_file_name.substr(dot_pos);
        result_file_name = result_file_name.substr(0, dot_pos);

        dot_pos = result_file_name.find_last_of(".");
        if (dot_pos != std::string::npos) {
            pre_last_extension = result_file_name.substr(dot_pos);
            result_file_name = result_file_name.substr(0, dot_pos);
        }
    }

    std::transform(last_extension.begin(), last_extension.end(), last_extension.begin(), ::tolower);
    if (last_extension == GZ_FILE_EXTENSION) {
        type = FileReader::GZIP;
    }
    else {
        type = FileReader::RAW;
    }

    extension = pre_last_extension + last_extension;
    raw_name = result_file_name;
    correction_info_name = raw_name + CORRECTION_FILE_EXTENSION; 
    output_name = raw_name + CORRECTED_FILE_EXTENSION + extension;
}



void ReadFileData::addPrefix(const std::string& prefix) {
    correction_info_name = prefix + DIRECTORY_SEPARATOR + correction_info_name;
    output_name = prefix + DIRECTORY_SEPARATOR + output_name;
}



std::string ReadFileData::getCorrectionInfoName(const std::size_t fileNumber) const {
    if (fileNumber == 0) {
        return output_name;
    }

    std::ostringstream error_correction_info_file_name_stream;
    error_correction_info_file_name_stream << correction_info_name;
    error_correction_info_file_name_stream << fileNumber;

    return error_correction_info_file_name_stream.str();
}



//----------------------------------------------------------------------
// Parses program arguments.
//----------------------------------------------------------------------

void C_arg::read_args() {
    std::cout << "RECKONER ver. " << VERSION << std::endl << std::endl;
    if (num_args == 1 || (num_args == 2 && (strcmp(args[1], "-help") == 0))) {
        print_usage();
        exit(EXIT_SUCCESS);
    }

    bool lastArgumentsGiven = false;

    // to verify, if user has not given arguments configuring long k-mers verification when such a verification is not enabled
    bool accept_given = false;
    bool longratio_given = false;

    bool kmc_memory_given = false;

    //--------------------------------------------------
    // parse arguments
    //--------------------------------------------------
    for (int it_arg = 1; it_arg < num_args; it_arg++) {
        if (strcmp(args[it_arg], "-read") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                read_files_data.push_back(ReadFileData(args[it_arg + 1]));
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The read file name is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-prefix") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                prefix = args[it_arg + 1];
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The prefix is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-kmerlength") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                try {
                    kmer_length = std::stol(args[it_arg + 1]);
                }
                catch (...) {
                    std::cerr << "ERROR: The k-mer length must be numeric" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                
                is_kmer_length_user_defined = true;
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The k-mer length is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-genome") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                try {
                    genome_size = std::stol(args[it_arg + 1]);
                }
                catch (...) {
                    std::cerr << "ERROR: The genome size must be numeric" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                is_genome_size_user_defined = true;
                it_arg++;
            }
            else {
                std::cerr << "ERROR: Genome size is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-extend") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                try {
                    extend = std::stol(args[it_arg + 1]);
                }
                catch (...) {
                    std::cerr << "ERROR: The max extension must be numeric" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The max extension is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-threads") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                try {
                    n_threads = std::stol(args[it_arg + 1]);
                }
                catch (...) {
                    std::cerr << "ERROR: The number of threads must be numeric" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The number of threads is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-noindels") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            correct_indels = false;
        }
        else if (strcmp(args[it_arg], "-fixlength") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            change_headers_length_tag = true;
        }
        else if (strcmp(args[it_arg], "-longkmer") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            use_long_kmer = true;
        }
        else if (strcmp(args[it_arg], "-accept") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            accept_filtered_with_long_kmers = true;
            accept_given = true;
        }
        else if (strcmp(args[it_arg], "-longratio") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                try {
                    long_kmer_ratio = std::stof(args[it_arg + 1]);
                }
                catch (...) {
                    std::cerr << "ERROR: The long to short k-mer ratio must be numeric" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The long to short k-mer ratio is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
            if (long_kmer_ratio <= 1.0) {
                std::cerr << "ERROR: Long k-mer ratio has to be > 1.0." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
            longratio_given = true;
        }
        else if (strcmp(args[it_arg], "-mark") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            mark_corrected = true;
        }
        else if (strcmp(args[it_arg], "-verbose") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            verbose = true;
        }
        // Options for debugging
        else if (strcmp(args[it_arg], "-kmcmemory") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            if (it_arg <= num_args - 2) {
                try {
                    kmc_memory = std::stol(args[it_arg + 1]);
                }
                catch (...) {
                    std::cerr << "ERROR: The limit of k-mer counting memory usage must be numeric" << std::endl << std::endl;
                    exit(EXIT_FAILURE);
                }
                kmc_memory_given = true;
                it_arg++;
            }
            else {
                std::cerr << "ERROR: The limit of k-mer counting memory usage is not specified." << std::endl << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(args[it_arg], "-reuse") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            reuse_kmc_db = true;
        }
        else if (strcmp(args[it_arg], "-kmcram") == 0) {
            if (lastArgumentsGiven) {
                params_at_the_beginning();
            }
            kmc_ram = true;
            kmc_memory = MAX_KMC_RAM_MEMORY_USAGE;
        }
        else if (args[it_arg][0] == '-') {
            std::cerr << "ERROR: Illegal option " << args[it_arg] << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }
        else {
            lastArgumentsGiven = true;
            read_files_data.push_back(ReadFileData(args[it_arg]));
        }
    }


    //--------------------------------------------------
    // check options
    //--------------------------------------------------
    if (read_files_data.empty()) {
        std::cerr << "ERROR: No read file name is specified." << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (prefix.empty()) {
        prefix = ".";
    }

    for (auto& it : read_files_data) {
        it.addPrefix(prefix);
    }

    if (is_kmer_length_user_defined && kmer_length < MIN_KMER_LENGTH) {
        std::cerr << "ERROR: k-mer length should be >= " << MIN_KMER_LENGTH << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (extend < 1) {
        std::cerr << "ERROR: The max extension should be >= 1" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (n_threads < 1) {
        std::cerr << "ERROR: Number of threads should be >= 1" << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!use_long_kmer) {
        if (accept_given) {
            std::cerr << "ERROR: -accept parameter can be given only when -longkmer is enabled." << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }
        if (longratio_given) {
            std::cerr << "ERROR: -longratio parameter can be given only when -longkmer is enabled." << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!correct_indels && change_headers_length_tag) {
        std::cerr << "ERROR: -fixlength parameter cannot be given with -noindels." << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (kmc_ram && kmc_memory_given) {
        std::cerr << "ERROR: -kmcmemory parameter cannot be given with -kmcram." << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    if (use_long_kmer && is_kmer_length_user_defined) {
        long_kmer_length = static_cast<std::size_t>(long_kmer_ratio * kmer_length);
    }

    // set auxiliary files names
    log_file_name = prefix + DIRECTORY_SEPARATOR + read_files_data.front().raw_name + LOG_FILE_EXTENSION;
    kmc_determine_params_database_name = prefix + DIRECTORY_SEPARATOR + read_files_data.front().raw_name + DETERMINE_PARAMS_EXTENSION;
    kmc_database_name = prefix + DIRECTORY_SEPARATOR + read_files_data.front().raw_name + TEMP_EXTENSION;
    kmc_filtered_database_name = prefix + DIRECTORY_SEPARATOR + read_files_data.front().raw_name;
    if (use_long_kmer) {
        kmc_long_database_name = prefix + DIRECTORY_SEPARATOR + read_files_data.front().raw_name + LONG_KMERS_EXTENSION;
    }
    kmc_list_file_name = prefix + DIRECTORY_SEPARATOR + read_files_data.front().raw_name + LIST_FILE_EXTENSION;

    //--------------------------------------------------
    // file names
    //--------------------------------------------------

    for (std::size_t it = 0; it < read_files_data.size(); ++it) {
        std::ifstream f_tmp;
        f_tmp.open(read_files_data[it].input_name);

        if (f_tmp.is_open() == false) {
            std::cerr << "ERROR: Cannot open " << read_files_data[it].input_name << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }

        f_tmp.close();
    }

    if (!RunExternal::createDirectory(prefix)) {
        std::cerr << "ERROR: Cannot create output directory (" << prefix << ")." << std::endl;
        exit(EXIT_FAILURE);
    }

    Log::open_log_file(log_file_name);
    Log::get_stream() << "RECKONER ver. " << VERSION << std::endl << std::endl;
    C_log c_log(std::cout);

    //--------------------------------------------------
    // print options
    //--------------------------------------------------
    c_log << "##################################################" << std::endl;
    c_log << "INPUT DATA SUMMARY" << std::endl;
    c_log << "##################################################" << std::endl;

    for (std::size_t it = 0; it < read_files_data.size(); ++it) {
        c_log << "Read file " << (it + 1) << std::endl;
        c_log << "      Name                    : " << read_files_data[it].input_name << std::endl;
        std::string compressed = "no";
        if (read_files_data[it].type == FileReader::GZIP) {
            compressed = "yes (GZIP)";
        }
        c_log << "      Compressed              : " << compressed << std::endl;
    }
    c_log << "Log file name                 : " << log_file_name << std::endl;
    if (is_kmer_length_user_defined) {
        c_log << "K-mer length                  : " << kmer_length << std::endl;
    }
    else {
        c_log << "K-mer length                  : " << "not defined (will be determined automatically)" << std::endl;
    }
    if (use_long_kmer) {
        if (is_kmer_length_user_defined) {
            c_log << "Long k-mer length             : " << long_kmer_length << std::endl;
        }
        else {
            c_log << "Long k-mer length             : " << "not defined (will be determined automatically)" << std::endl;
        }
    }
    c_log << "Number of threads             : " << n_threads << std::endl;
}



//----------------------------------------------------------------------
// Prints arguments usage.
//----------------------------------------------------------------------

void C_arg::print_usage() {
    std::cout << std::endl;
    std::cout << "USAGE: " << args[0] << " {ARGUMENTS} {FASTQ files}" << std::endl;
    std::cout << std::endl;

    std::cout << "Main arguments (sufficient in a general use):" << std::endl;
    std::cout << "-help                   - print this text" << std::endl;
    std::cout << "-read <FASTQ file>      - input FASTQ(.gz) read file name, can be passed many times" << std::endl;
    std::cout << "-kmerlength <K>         - k-mer length" << std::endl;

    std::cout << "Supporting:" << std::endl;
    std::cout << "-genome <G>             - approximate genome size, not mandatory: when -kmerlength is not given, it can help k-mer length determination" << std::endl;

    std::cout << "Quality tunning:" << std::endl;
    std::cout << "-longkmer               - verify corrections with long k-mers (worthful for de novo assembly, requires more RAM and CPU time), it enables:" << std::endl;
    std::cout << "    -accept             - if it is not possible to verify some read, keep the best of its possible corrections (default: do not correct the read)" << std::endl;
    std::cout << "    -longratio <RATIO>  - long to normal k-mer length ratio, must be > 1.0, default " << DEFAULT_LONG_KMER_RATIO << std::endl;
    std::cout << "-noindels               - do not correct insertions and deletions (substitutions correction only)" << std::endl;
    std::cout << "-extend <N>             - max extend length, default " << MAX_EXTENSION << std::endl;

    std::cout << "Formal issues:" << std::endl;
    std::cout << "-prefix <DIRECTORY>     - output directory, default current directory (.)" << std::endl;
    std::cout << "-threads <N>            - number of correcting threads, default number of available virtual cores" << std::endl;
    std::cout << "-mark                   - mark corrected reads by adding an information to their headers" << std::endl;
    std::cout << "-fixlength              - if the read headers contain read length (e.g. \"length=100\"), fix this value after indel correction (cannot be used with -noindels)" << std::endl;
    
    std::cout << "Useful for debugging:" << std::endl;
    std::cout << "-kmcmemory <N>          - k-mer counting memory consumption limit in GB, default " << DEFAULT_MAX_KMC_MEMORY_USAGE << std::endl;
    std::cout << "-kmcram                 - count k-mers in RAM-only mode" << std::endl;
    std::cout << "-reuse                  - do not remove a KMC database; if a proper database exists; if proper database exists, use it" << std::endl;
    std::cout << "-verbose                - return more statistics" << std::endl;
    

    std::cout << std::endl << "Instead of using -read option the input FASTQ files can be specified" << std::endl;
    std::cout << "at the end of the command line." << std::endl;

    std::cout << std::endl << "K-mer length, if not given, will be determined automatically. Delivering genome size" << std::endl;
    std::cout << "can facilitate and accelerate this determination." << std::endl;

    std::cout << std::endl;
}



//----------------------------------------------------------------------
// Handles wrong read files placement.
//----------------------------------------------------------------------

void C_arg::params_at_the_beginning() {
    std::cerr << "ERROR: Input files list should be specified at the end of the command line." << std::endl;
    std::cerr << "       Alternatively use -read option." << std::endl;
    exit(EXIT_FAILURE);
}
