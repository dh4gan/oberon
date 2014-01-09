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
	Planet(string &namestring, string &typestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel);
	Planet(string &namestring, string &typestring, double &m, double &rad, double semimaj, double ecc, double inc,
			double longascend, double argper, double meananom, double G, double totalMass);
	virtual ~Planet();

	/* Accessors */

	// Standard cloning method
	virtual Planet* Clone() { return new Planet(*this); }


};

#endif /* PLANET_H */

