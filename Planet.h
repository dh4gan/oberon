/*
 * Planet.h
 *
 *  Created on: Nov 8, 2012
 *      Author: dh4gan
 */

#include <iostream>
#include "Body.h"

using namespace std;

#ifndef PLANET_H_
#define PLANET_H_

class Planet: public Body {
public:

	Planet();
	Planet(string &namestring, string &typestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &alb);
	Planet(string &namestring, string &typestring, double &m, double &rad, double semimaj, double ecc, double inc,
			double longascend, double argper, double meananom, double G, double totalMass, double &alb);
	virtual ~Planet();

	/* Accessors */

	// Standard cloning method
	virtual Planet* Clone() { return new Planet(*this); }

	void setEquilibriumTemperature(double temp){temperature=temp;}
	double getEquilibriumTemperature(){return temperature;}

	void setReflectiveLuminosity(double lum){reflectiveLuminosity = lum;}
	void setAlbedo (double alb){albedo = alb;}
	double getAlbedo(){return albedo;}

	void calcLuminosity();

	void setLuminosity(double lum){luminosity = lum;}
	double getLuminosity(){return luminosity;}

protected:
	double temperature;
	double albedo;
	double reflectiveLuminosity;
	double luminosity;
};

#endif /* PLANET_H */

