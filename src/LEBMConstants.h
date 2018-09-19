/*
 * LEBMConstants.h
 *
 *  Created on: Sep 19, 2018
 *      Author: dhf
 *      Constants relevant to the 1D Latitudinal Energy Balance Model (LEBM)
 */

#ifndef LEBMCONSTANTS_H_
#define LEBMCONSTANTS_H_

const double freeze = 273.0; // Freezing point of water
const double boil = 373.0; // Boiling point of water

// Minimum and maximum temperature ranges for humans (without shelter)
const double humanTooCold = 277.0
const double humanTooHot = 323.0

// Surface heat capacity coefficients
const double C_land = 5.25e9; //differs by factor 1000 to W&K due to units...
const double C_ocean = 40.0 * C_land;


#endif /* LEBMCONSTANTS_H_ */
