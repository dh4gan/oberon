/*
 * Star.h
 *
 *  Created on: Nov 8, 2012
 *      Author: dh4gan
 */

#include "Body.h"

using namespace std;

#ifndef STAR_H_
#define STAR_H_

class Star: public Body {
public:

	Star();
	Star(string &namestring, string &typestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &lum);
	Star(string &namestring, string &typestring, double &m, double &rad, double semimaj, double ecc, double inc,
			double longascend, double argper, double meananom, double G, double totalMass, double &lum);
	virtual ~Star();

	void setLuminosity(double lum){luminosity = lum;}
	double getLuminosity() {return luminosity;}

	void calcMainSequenceLuminosity(){luminosity = pow(mass,4);}

	// Standard cloning method
	virtual Star* Clone() { return new Star(*this); }

protected:

	double luminosity; // Luminosity of Star IN SOLAR UNITS

};

#endif /* STAR_H */

