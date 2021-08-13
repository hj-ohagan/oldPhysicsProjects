//
//  ising.h
//  ising
//
//  Created by Henry O'Hagan on 01/12/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

#ifndef ising_ising_h
#define ising_ising_h

#include "macrostate.h"
#include <ctime>

void simulate_constantbfield(std::ofstream& ,double , std::ofstream& );
void simulate_constanttemp(std::ofstream&,std::ofstream&  ,double );
void printtime(time_t);
void addresults(int , std::vector<std::vector<double>>& , macrostate, bool);
void printresults(int ,int ,int ,std::vector<std::vector<double>>  , macrostate  ,std::ofstream& , bool);
void test(std::ofstream&);
#endif
