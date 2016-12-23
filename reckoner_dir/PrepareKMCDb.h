/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* 0.2.1
*
*/

#ifndef PREPAREKMCDB_H
#define	PREPAREKMCDB_H



#include "parse_args.hpp"


class PrepareKMCDb {
public:
    void run(const C_arg& c_inst_args);
};

#endif	/* PREPAREKMCDB_H */
