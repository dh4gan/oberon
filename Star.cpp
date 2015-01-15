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
    luminosity = 0.0;
    calcMainSequenceLuminosity();
}
Star::Star(string &namestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &lum) :
	Body(namestring, m, rad, pos, vel) {

	type = "Star";
	if(lum > 0.0)
	    {
	    luminosity = lum;
	    }
	else
	    {
	    calcMainSequenceLuminosity();
	    }
}

Star::Star(string &namestring, double &m, double &rad, double semimaj, double ecc, double inc,
			double longascend, double argper, double meananom, double G, double totalMass, double &lum):
			Body(namestring, m, rad, semimaj,ecc,inc,longascend,argper,meananom,G,totalMass) {
	type = "Star";
	if(lum > 0.0)
    	    {
    	    luminosity = lum;
    	    }
    	else
    	    {
    	    calcMainSequenceLuminosity();
    	    }

}


Star::~Star() {
}

