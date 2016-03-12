/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#include "time.hpp"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "PerformTesting.h"
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

    // corrent errors in reads
    performTesting.perform_correction(c_inst_args, c_inst_time, check_inputs);

    return (EXIT_SUCCESS);
}
