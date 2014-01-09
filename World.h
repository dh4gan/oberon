/*
 * World.h
 *
 *  Created on: Jan 9, 2014
 *      Author: dh4gan
 */

#include <iostream>
#include "Body.h"
#include "Star.h"

using namespace std;

#ifndef WORLD_H_
#define WORLD_H_

class World: public Body
    {
public:

    World();
    World(string &namestring, string &typestring, double &m, double &rad,
	    Vector3D &pos, Vector3D &vel, int &n, double &obliq, double &winter,
	    double &ocean, double &T0);
    World(string &namestring, string &typestring, double &m, double &rad,
	    double semimaj, double ecc, double inc, double longascend,
	    double argper, double meananom, double G, double totalMass,  int &n, double &obliq, double &winter,
	    double &ocean, double &T0);
    virtual ~World();

    /* Accessors */

    void setInsolationZero() {insol.assign(nPoints1,0.0);}

    // Calculation Methods

    void initialiseLEBM();
    void calcInsolation(Star star); // TODO
    void checkForEclipses(Star star, vector<Body*> bodies); // TODO

    void calcAlbedo();
    void calcHeatCapacity();
    void calcIce();
    void calcOpticalDepth();
    void calcCooling();
    void calcNetHeating();
    void calcHabitability(double &minT, double &maxT);

    void integrate(double &dt); //TODO

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
    double landFraction;
    double initialTemperature;
    double diffusion;

    int nPoints,nPoints1;

    double nFloat;

    vector<double> lat, x, deltax;
    vector<double> T, T_old, tau;
    vector<double> iceFraction, C, hab;
    vector<double> infrared, Q, albedo, insol;

    };

#endif /* WORLD_H */

