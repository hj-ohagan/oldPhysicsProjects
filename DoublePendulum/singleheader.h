//
//  singleheader.h
//  CompPhysSinglePendulum
//
//  Created by Henry O'Hagan on 04/11/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

#ifndef CompPhysSinglePendulum_singleheader_h
#define CompPhysSinglePendulum_singleheader_h

/*-----------------------------------------------------CLASSES-------------------------------------------------------*/

#include <fstream>
#include <cmath>
#include "matrix.h"
#include <complex>


/*-------------------------------------------------FUNCTION_HEADERS--------------------------------------------------*/

matrix euler_updater_initialise(double , double );
matrix leapfrog_updater_initialise(double , double );
matrix rk_updater_initialise(double, double);

vect analytic(double, double , double);

std::vector<std::complex<double>> deigen(double, double);
std::vector<std::complex<double>> eeigen(double, double);
std::vector<std::complex<double>> rkeigen(double, double);

double hdifference(vect , vect );
double energy(vect );

void single_step_updater( double , double , vect , std::ofstream&, char );
void leapfrog(double , double , vect , std::ofstream& );
void criticalstep(std::ofstream&, char);

#endif /* defined(__SinglePendulum__File__) */

