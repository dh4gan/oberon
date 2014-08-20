/*
 * Constants.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dhf
 *      NOTE ALL CONSTANTS ARE IN SI UNITS
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

double pi = 3.141592;
double Gsi = 6.67e-11;
double Gmau = 1.0;  // Value of G for solar mass-AU units (time units...2pi units= 1 year)
double AU = 1.496e11;
double msol = 1.99e30;
double lsol = 3.826e26;
double rsol = 6.955e8;

double fluxsol = lsol/(4.0*pi*AU*AU);
double fluxsolcgs = fluxsol*1000;
//double fluxsol = 1;



#endif /* CONSTANTS_H_ */
