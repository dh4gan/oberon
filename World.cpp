/*
 * World.cpp
 *
 *  Created on: Jan 9 2014
 *      Author: dh4gan
 */

#include "World.h"

World::World() :
	Body()
    {
    }
World::World(string &namestring, string &typestring, double &m, double &rad,
	Vector3D &pos, Vector3D &vel, int &n, double &obliq, double &winter,
	    double &ocean, double &T0) :
	Body(namestring, typestring, m, rad, pos, vel)
    {
    nPoints = n;
    obliquity = obliq;
    winterSolstice = winter;
    oceanFraction = ocean;
    initialTemperature = T0;



    }
World::World(string &namestring, string &typestring, double &m, double &rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, int &n, double &obliq, double &winter,
	    double &ocean, double &T0) :
	Body(namestring, typestring, m, rad, semimaj, ecc, inc, longascend,
		argper, meananom, G, totalMass)
    {
    nPoints = n;
    obliquity = obliq;
    winterSolstice = winter;
    oceanFraction = ocean;
    initialTemperature = T0;

    }
World::~World()
    {
    }

