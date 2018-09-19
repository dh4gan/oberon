/*
 * Constants.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dhf
 *      Defines standard physical constants
 *      NOTE ALL CONSTANTS ARE IN SI UNITS UNLESS OTHERWISE STATED
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double pi = 3.141592653;
const double twopi = 2.0*pi;
const double piby2 = 0.5*pi;

const double year = 3.1556926e7; // year in seconds
const double degToRad = pi/180.0;
const double radToDeg = 1.0/degToRad;

const double AU = 1.496e11;
const double msol = 1.99e30;
const double lsol = 3.8626e26;
const double rsol = 6.955e8;

const double solradToAU = 0.00467580213904;

const double Gsi = 6.67e-11;
const double Gmau = 1.0;  // Value of G for solar mass-AU units (time units...2pi units= 1 year)
const double Gmau_day = 2.959e-4; // Value of G for (Msol, AU, day) unit system)

const double fluxsol = lsol/(4.0*pi*AU*AU); // Solar constant
const double fluxsolcgs = fluxsol*1000; // Solar constant in cgs units

const double hplanck = 6.626e-34;
const double c = 2.9e8;
const double k_B = 1.38e-23;
const double sigma_SB = 5.67e-5; //erg s^-1 cm^-2 K^-4
const double c_cgs = 3.0e10;
const double c_mau = c_cgs*year/(twopi*AU);

#endif /* CONSTANTS_H_ */
