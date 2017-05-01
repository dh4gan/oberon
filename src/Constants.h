/*
 * Constants.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dhf
 *      NOTE ALL CONSTANTS ARE IN SI UNITS
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double pi = 3.141592;
const double twopi = 2.0*pi;
const double piby2 = 0.5*pi;

const double year = 3.15e7; // year in seconds
const double degToRad = pi/180.0;
const double radToDeg = 1.0/degToRad;
double Gsi = 6.67e-11;
double Gmau = 1.0;  // Value of G for solar mass-AU units (time units...2pi units= 1 year)
double Gmau_day = 2.959e-4; // Value of G for (Msol, AU, day) unit system)
double AU = 1.496e11;
double msol = 1.99e30;
double lsol = 3.826e26;
double rsol = 6.955e8;

double fluxsol = lsol/(4.0*pi*AU*AU);
double fluxsolcgs = fluxsol*1000;
//double fluxsol = 1;

const double c_cgs = 3.0e10;
const double c_mau = c_cgs*year/(twopi*AU);

const double CO2Earth = 3.3e-4;
const double kactive = 0.09;
const double krun = 0.045;

#endif /* CONSTANTS_H_ */
