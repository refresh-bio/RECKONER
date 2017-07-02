/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.1
 * 
 */

#include "time.hpp"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "PerformTasks.h"
#include <vector>



//----------------------------------------------------------------------
// Main function.
//----------------------------------------------------------------------

int main(int argc, char** argv) {
    const C_arg c_inst_args(argc, argv);

    PerformTasks performTasks(c_inst_args);
    performTasks.perform_check_read();
    performTasks.perform_parameter_determination();
    performTasks.perform_kmcdb_creation();
    performTasks.perform_correction();
    performTasks.finalize();

    return (EXIT_SUCCESS);
}
