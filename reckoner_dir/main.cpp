/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.2
 * 
 */

#include "DetermineParameters.hpp"
#include "time.hpp"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "PerformTesting.h"
#include "PrepareKMCDb.h"
#include <vector>



//----------------------------------------------------------------------
// Main function.
//----------------------------------------------------------------------

int main(int argc, char** argv) {
    // construct initial classes
    C_time c_inst_time;
    C_arg c_inst_args(argc, argv, c_inst_time);

    std::vector<C_check_read> check_inputs;

    PerformTasks performTesting;

    // check input read files
    performTesting.perform_check_read(c_inst_args, c_inst_time, check_inputs);

    DetermineParameters determineParameters(c_inst_args, c_inst_time, check_inputs);
    determineParameters.perform_determine_parameters();

    PrepareKMCDb prepareKMCDb;
    prepareKMCDb.run(c_inst_args);

    std::cout << std::endl;
    std::cout << "##################################################" << std::endl;
    std::cout << "CORRECTING ERRORS" << std::endl;
    std::cout << "##################################################" << std::endl;
    Log::get_stream() << std::endl;
    Log::get_stream() << "##################################################" << std::endl;
    Log::get_stream() << "CORRECTING ERRORS" << std::endl;
    Log::get_stream() << "##################################################" << std::endl;

    // correct errors in reads
    performTesting.perform_correction(c_inst_args, c_inst_time, check_inputs);

    return (EXIT_SUCCESS);
}
