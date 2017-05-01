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
	Star(string &namestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel, double &lum);
	Star(string &namestring, double &m, double &rad, double semimaj, double ecc, double inc,
			double longascend, double argper, double meananom, double G, double totalMass, double &lum);
	virtual ~Star();

	void setLuminosity(double lum){luminosity = lum;}
	double getLuminosity() {return luminosity;}

	double getfVisible(){return fVisible;}
	double getfIR(){return fIR;}

	void calcMainSequenceLuminosity(){luminosity = pow(mass,4);}

	// Standard cloning method
	virtual Star* Clone() { return new Star(*this); }

protected:

	double luminosity; // Luminosity of Star IN SOLAR UNITS
	double fVisible; // Fraction of flux in visible
	double fIR; // Fraction of flux in IR

};

#endif /* STAR_H */

