/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 2.0
*
*/

#ifndef PREPAREKMCDB_H
#define PREPAREKMCDB_H

#include "parse_args.hpp"
#include "time.hpp"
#include "Log.h"



class PrepareKMCDb {
private:
    const C_arg& c_inst_args;
    C_log c_err;

public:
    PrepareKMCDb(const C_arg& _c_inst_args) : c_inst_args(_c_inst_args), c_err(std::cerr) {}

    void countKmers();
    void countLongKmers();
    unsigned determineCutoff();
    void filterKmers(unsigned cutoff);

    void removeKMCDatabase();
    void removeFilteredKMCDatabase();
};

#endif /* PREPAREKMCDB_H */
