/*
 * RECKONER - Read Error Corrector Based on KMC
 * 
 * This software is distributed under GNU GPL 3 license.
 * 
 * Authors: Yun Heo, Maciej Dlugosz
 * Version: 0.1
 * 
 */

#ifndef PERFORMTESTING_H
#define	PERFORMTESTING_H

#include "parse_args.hpp"
#include "check_inputs.hpp"
#include <vector>
#include <string>

class PerformTasks {
public:
    void perform_check_read(const C_arg& c_inst_args, C_time& c_inst_time, std::vector<C_check_read>& check_inputs);
    void summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time);
    void perform_correction(const C_arg& c_inst_args, C_time& c_inst_time, std::vector<C_check_read>& check_inputs);
private:
    std::string remove_new_line(std::string in_string);
};

#endif	/* PERFORMTESTING_H */

