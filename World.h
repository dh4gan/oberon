/*
 * World.h
 *
 *  Created on: Jan 9, 2014
 *      Author: dh4gan
 */

#include <iostream>
#include "Body.h"

using namespace std;

#ifndef WORLD_H_
#define WORLD_H_

class World: public Body
    {
public:

    World();
    World(string &namestring, string &typestring, double &m, double &rad,
	    Vector3D &pos, Vector3D &vel, int &n, double &obliquity, double &winter,
	    double &ocean, double &T0);
    World(string &namestring, string &typestring, double &m, double &rad,
	    double semimaj, double ecc, double inc, double longascend,
	    double argper, double meananom, double G, double totalMass,  int &n, double &obliquity, double &winter,
	    double &ocean, double &T0);
    virtual ~World();

    /* Accessors */


    void initialiseLEBM();

    // Standard cloning method
    virtual World* Clone()
	{
	return new World(*this);
	}

protected:
    double rotationPeriod;
    double obliquity;
    double winterSolstice;
    double oceanFraction;
    double initialTemperature;

    int nPoints;

    vector<double> lat, x, deltax;
    vector<double> T, T_old;
    vector<double> f_ice, hab;
    vector<double> infrared, Q, albedo, insol;

    };

#endif /* WORLD_H */

