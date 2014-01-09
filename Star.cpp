/*
 * Star.cpp
 *
 *  Created on: 13 Sep 2012
 *      Author: dhf
 */

#include "Star.h"
#include "RadiationConstants.h"
#include <iostream>
/* Constructors and Destructor */

Star::Star() :
	Body() {
}
Star::Star(string &namestring, string &typestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &lum) :
	Body(namestring, typestring, m, rad, pos, vel) {
	luminosity = lum;
}

Star::Star(string &namestring, string &typestring, double &m, double &rad, double semimaj, double ecc, double inc,
			double longascend, double argper, double meananom, double G, double totalMass, double &lum):
			Body(namestring, typestring, m, rad, semimaj,ecc,inc,longascend,argper,meananom,G,totalMass) {
    luminosity = lum;
}


Star::~Star() {
}

