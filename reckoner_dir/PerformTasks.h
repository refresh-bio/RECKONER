/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 1.1
 * 
 */

#ifndef PERFORMTESTING_H
#define	PERFORMTESTING_H

#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "Log.h"
#include <vector>
#include <string>



class PerformTasks {
public:
    PerformTasks(const C_arg& _c_inst_args) : c_inst_args(_c_inst_args), f_log(Log::get_stream()) {};

    void perform_check_read();
    void perform_parameter_determination();
    void perform_kmcdb_creation();
    void perform_correction();
    void finalize();

private:
    C_time c_inst_time;
    const C_arg& c_inst_args;
    std::vector<C_check_read> check_inputs;

    std::ofstream& f_log;

    void summarize_outputs();
};

#endif	/* PERFORMTESTING_H */

