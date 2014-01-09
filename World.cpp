/*
 * World.cpp
 *
 *  Created on: Jan 9 2014
 *      Author: dh4gan
 */

#include "World.h"

double pi = 3.141592;
double freeze = 273.0;
double boil = 373.0;
double sigma_SB = 5.67e-5;

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
    landFraction = 1.0-oceanFraction;
    initialTemperature = T0;
    nFloat = float(nPoints);

    nPoints1 = nPoints+1;

    initialiseLEBM();

    }
World::World(string &namestring, string &typestring, double &m, double &rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, int &n,
	double &obliq, double &winter, double &ocean, double &T0) :
	Body(namestring, typestring, m, rad, semimaj, ecc, inc, longascend,
		argper, meananom, G, totalMass)
    {
    nPoints = n;
    obliquity = obliq;
    winterSolstice = winter;
    oceanFraction = ocean;
    landFraction = 1.0-oceanFraction;
    initialTemperature = T0;
    nFloat = float(nPoints);

    nPoints1 = nPoints+1;

    initialiseLEBM();

    }
World::~World()
    {
    }


void World::initialiseLEBM()
    {
    /*
     * Written 9/1/14 by dh4gan
     * Sets up the LEBM vectors for computation
     */

    // Set length of vectors

   lat.resize(nPoints1,0.0);
   x.resize(nPoints1,0.0);
   deltax.resize(nPoints1,0.0);
   iceFraction.resize(nPoints1,0.0);
   hab.resize(nPoints1,0.0);
   infrared.resize(nPoints1,0.0);
   Q.resize(nPoints1,0.0);
   albedo.resize(nPoints1,0.0);
   insol.resize(nPoints1,0.0);
   tau.resize(nPoints1,0.0);
   C.resize(nPoints1,0.0);

   // Set temperature equal to initial temperature
   T.resize(nPoints1,initialTemperature);
   T_old.resize(nPoints1, 0.0);

   // Set up latitudes etc

   double dlat = pi/nFloat;

   for (int i=0; i< nPoints1; i++)
       {
       lat[i] = -pi/2.0 + i*dlat;
       x[i] = sin(lat[i]);
       deltax[i] = cos(lat[i])*dlat;
       }

   // Set up diffusion constant TODO

    }


void World::calcAlbedo()
    {
    /*
     * Written 9/1/14 by dh4gan
     * Simple calculation of the albedo as a function of latitude according to Temperature
     *
     */


    for (int i=0; i< nPoints1; i++)
	{
	albedo[i] = 0.525 - 0.245* tanh((T[i]-freeze+5.0)/5.0);
	}

    }

void World::calcHeatCapacity()
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculation of Heat Capacity of the atmosphere as a function of the local temperature
     *
     */
    double C_ice;
    double C_land = 5.25e9;
    double C_ocean = 40.0 * C_land;

    for (int i = 0; i < nPoints1; i++)
	{
	if (T[i] >= freeze)
	    {
	    C_ice = 0.0;
	    }
	else if (T[i] < freeze and T[i] > freeze - 10.0)
	    {
	    C_ice = 9.2 * C_land;
	    }

	else if (T[i] <= freeze - 10)
	    {
	    C_ice = 2.0 * C_land;
	    }

	C[i] = landFraction * C_land
		+ oceanFraction*(iceFraction[i] * C_ice
		+ (1.0 - iceFraction[i]) * C_ocean);

	}

    }

void World:: calcIce()
    {

    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Ice Fraction given the temperature T
     *
     */

    for (int i=0; i< nPoints1; i++)
	{
	iceFraction[i] = 1.0 - exp(-(freeze-T[i])/10.0);
	}

    }

void World::calcOpticalDepth()
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the optical depth of the atmosphere
     */

    for (int i=0; i<nPoints1; i++)
	{
	tau[i] = 0.79*pow(T[i]/freeze,3);
	}
    }

void World::calcCooling()
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Infrared Cooling as a function of optical depth and temperature
     *
     */

    for(int i=0; i<nPoints1; i++)
	{
	infrared[i] = sigma_SB*pow(T[i],4)/(1.0+0.75*tau[i]);
	}

    }

void World::calcNetHeating()
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Net Heating in the LEBM system
     *
     */

    for (int i=0; i< nPoints1; i++)
	{
	Q[i] = insol[i]*(1.0-albedo[i]) - infrared[i];
	}
    }

void World::calcHabitability(double &minT, double &maxT)
    {
    /*
     * Written 9/1/14 by dh4gan
     *
     * Calculates a habitability index based on a specified set of
     * minimum and maximum temperatures for life
     * for life
     */

    for (int i=0; i< nPoints1; i++)
	{
	if(T[i]>=minT and T[i]<=maxT) hab[i]=1;
	}
    }
