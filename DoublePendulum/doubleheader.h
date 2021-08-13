//
//  doubleheader.h
//  doublependulum
//
//  Created by Henry O'Hagan on 06/11/2013.
//  Copyright (c) 2013 Henry O'Hagan. All rights reserved.
//

#ifndef doublependulum_doubleheader_h
#define doublependulum_doubleheader_h

/*-----------------------------------------------------CLASSES-------------------------------------------------------*/
#include "/Users/henryohagan/Documents/Compphys/ProjectA/1/CompPhysSinglePendulum/CompPhysSinglePendulum/matrix.h"
#include <fstream>
#include <cmath>
#include <complex>

/*-------------------------------------------------FUNCTION_HEADERS--------------------------------------------------*/
matrix rk_updater_initialise(double, double, double);

void rungekutta(  double , vect , std::ofstream&, std::ofstream& );

double energy(vect);
double hdifference(vect, vect);


#endif
