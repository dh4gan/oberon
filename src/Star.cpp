/*
 * Star.cpp
 *
 *  Created on: 13 Sep 2012
 *      Author: dhf
 */

#include "Star.h"
#include "Constants.h"
#include "CSCycleConstants.h"
#include <iostream>
#include <fstream>
#include <sstream>

/* Constructors and Destructor */

Star::Star() :
	Body() {
    luminosity = 0.0;

    calcMainSequenceLuminosity();
}
Star::Star(string &namestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &lum, string &spec) :
	Body(namestring, m, rad, pos, vel) {

	type = "Star";

	spectralType = spec;
	loadAlbedoCoefficients();
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
			double longascend, double argper, double meananom, double G, double totalMass, double &lum, string &spec):
			Body(namestring, m, rad, semimaj,ecc,inc,longascend,argper,meananom,G,totalMass) {
	type = "Star";

                
	spectralType = spec;

	loadAlbedoCoefficients();

	if(lum > 0.0)
    	    {
    	    luminosity = lum;
    	    }
    	else
    	    {
    	    calcMainSequenceLuminosity();
    	    }

}


Star::Star(parFile &input, int &bodyIndex, double &G):
Body(input,bodyIndex,G) {
    type = "Star";
    printf("Spectral Type read as %s \n",spectralType.c_str());
    if(input.getBoolVariable("CarbonateSilicateCycle"))
    {
    spectralType = input.getStringVariable("SpectralType",bodyIndex);
        printf("Spectral Type read as %s \n",spectralType.c_str());
    loadAlbedoCoefficients();
    }
    if(input.getDoubleVariable("Luminosity",bodyIndex) > 0.0)
    {
        luminosity = input.getDoubleVariable("Luminosity",bodyIndex);
    }
    else
    {
        calcMainSequenceLuminosity();
    }
    
}

Star::~Star() {
}

void :: Star::loadAlbedoCoefficients()
    /*
     * Written 21/07/2017 by dh4gan
     * Loads the albedo fitting function coefficients from values stored in Constants.h
     * (also loads the visibility fractions)
     */


    {

    if (spectralType == "F" or spectralType == "f")
	{
	fVisible = fStarFVisible;
	hotAlbedoCoefficients = fStarHotAlbedo;
	coldAlbedoCoefficients = fStarColdAlbedo;
	}

    else if (spectralType == "G" or spectralType == "g")

	{
	fVisible = gStarFVisible;
	hotAlbedoCoefficients = gStarHotAlbedo;
	coldAlbedoCoefficients = gStarColdAlbedo;
	}
    else if (spectralType == "K" or spectralType == "k")
	{
	fVisible = kStarFVisible;
	hotAlbedoCoefficients = kStarHotAlbedo;
	coldAlbedoCoefficients = kStarColdAlbedo;
	}

    else if (spectralType == "M" or spectralType == "m")
	{
	fVisible = mStarFVisible;
	hotAlbedoCoefficients = mStarHotAlbedo;
	coldAlbedoCoefficients = mStarColdAlbedo;
	}


    fIR = 1.0 - fVisible;

    }


vector<double> Star::getAlbedoCoefficients(double temperature)

	/*
	 * Written 01/05/2017 by dh4gan
	 * Returns the albedo coefficients appropriate to the star's spectral
	 * type and the planet's surface temperature
	 *
	 */
	    {

	    if(temperature < 250.0)
		{
		return coldAlbedoCoefficients;

		}

	    else
		{
		return hotAlbedoCoefficients;
		}

	    }
