/*
 * Planet.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "Planet.h"

Planet::Planet() :
	Body()
    {
    }
Planet::Planet(string &namestring, string &typestring, double &m, double &rad,
	Vector3D &pos, Vector3D &vel) :
	Body(namestring, typestring, m, rad, pos, vel)
    {
    }
Planet::Planet(string &namestring, string &typestring, double &m, double &rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass) :
	Body(namestring, typestring, m, rad, semimaj, ecc, inc, longascend,
		argper, meananom, G, totalMass)
    {
    }
Planet::~Planet()
    {
    }

