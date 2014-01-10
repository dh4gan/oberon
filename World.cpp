/*
 * World.cpp
 *
 *  Created on: Jan 9 2014
 *      Author: dh4gan
 */

#include "World.h"
#include "Constants.h"

double freeze = 273.0;
double boil = 373.0;
double sigma_SB = 5.67e-5;

double q0 = lsol/(4.0*pi*AU*AU);

World::World() :
	Body()
    {
    }
World::World(string &namestring, string &typestring, double &m, double &rad,
	Vector3D &pos, Vector3D &vel, int n, double obliq, double winter,
	double ocean, double T0) :
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

   // Set up diffusion constant

   diffusion = 5.394e2 * rotationPeriod*rotationPeriod;

    }


void World::updateLEBM(vector<Body*> bodies, vector<double>eclipsefrac)
    {
    /*
     * Written 10/1/14 by dh4gan
     * This method is an umbrella method, using the other methods
     * written here to advance the LEBM model by dt, for convenient use by other methods
     *
     */
    setInsolationZero();

    for (int b=0; b< bodies.size(); b++)
	{
	if (bodies[b]->getType()=="Star")
	    {
	    calcInsolation(bodies[b],eclipsefrac[b]);
	    }
	}
    calcIce();
    calcHeatCapacity();
    calcOpticalDepth();
    calcAlbedo();
    calcCooling();
    calcNetHeating();
    calcHabitability(freeze,boil);

    calcLEBMTimestep();

    integrate();

    }

void World::updateLEBM(vector<Body*> bodies, vector<double> eclipsefrac, double dt)
    {

    dtLEBM = dt;
    updateLEBM(bodies,eclipsefrac);

    }


void World::calcInsolation(Body* star, double &eclipsefrac)
    {

    /*
     * Written by dh4gan, 10/1/14
     * Calculates the flux received at each latitude from input star
     * Given that eclipsefrac of the stellar light is blocked by other objects
     */

    vector<double> cos_H(nPoints1,0.0);

    // Calculate stellar declination

    // Get position of world relative to star
    Vector3D worldpos = this->getPosition();
    Vector3D starpos = star->getPosition();
    Vector3D pos = (starpos).relativeVector(worldpos);
    double magpos = pos.magVector();
    Vector3D unitpos = pos.unitVector();

    // Declination of the Sun - angle between planet's position vector and equator (at noon)

    Vector3D decVector(unitpos.elements[0], unitpos.elements[1],
			unitpos.elements[2]);

    // Rotate this vector if world has non-zero obliquity
    if (obliquity != 0.0) {
	decVector.rotateX(obliquity);
	}

    // Obtain declination angle

    double rdotn = unitpos.dotProduct(decVector);
    double declination = safeAcos(rdotn);

    double sind = sin(declination);
    double cosd = cos(declination);
    double tand = tan(declination);

    // Insolation prefactor depends on luminosity and separation only

    double lstar = star->getLuminosity();

    // Rotate in y-axis if world's winter solstice longitude non-zero TODO

    for (int i=0; i<nPoints1; i++)
	{

	// calculate the diurnally averaged hour angle

	cos_H[i] = -tan(lat[i])*tand;

	if(fabs(cos_H[i]) >1.0) cos_H[i] = cos_H[i]/fabs(cos_H[i]);

	double H = acos(cos_H[i]);

	insol[i] = insol[i]+q0*lstar*eclipsefrac/(pi*magpos*magpos)*(H*x[i]*sind + cos(lat[i])*cosd*sin(H));
	}


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

void World::calcLEBMTimestep()
    {
/*
 * Written 10/1/14 by dh4gan
 * Calculates the minimum timestep for the LEBM (in seconds)
 *
 */
    double Dplus, timestep;

    dtLEBM = 1.0e30;
    for (int i=0; i<nPoints1; i++)
	{
	if(i==nPoints)
	    {
	    Dplus = diffusion*(1.0-x[i]*x[i]);
	    }
	else
	    {
	    Dplus = diffusion*0.5*((1.0-x[i]*x[i]) +(1.0-x[i+1]*x[i+1]));
	    }

	if(i>0 and i<nPoints)
	    {
	    timestep = deltax[i]*deltax[i]*C[i]/(2.0*Dplus);
	    if (timestep < dtLEBM) {dtLEBM = timestep;}
	    }
	}


    if(dtLEBM==1.0e30)
	{
	cout << "ERROR in LEBM Timestep calculation " << endl;
	dtLEBM = -1.0;
	}

    }

void World::integrate()
    {
    /*
     * Written 10/1/14 by dh4gan
     * Integrates the diffusion equation to drive the system forward
     *
     */
    double Tminus1,Tplus1, Dplus,Dminus;
    double T1,dx,dx1,Fj;


    T_old = T;

    for(int i=0; i<nPoints1; i++)
	{
	if(i==0){
	    Tminus1 = T_old[i];
	    Dminus = diffusion*(1.0-x[i]*x[i]);
	}
	if(i/=1){
	    Tminus1 = T_old[i-1];
	    Dminus = 0.5*diffusion*((1.0-x[i-1]*x[i-1]) + (1.0-x[i]*x[i]));
	}

	if(i==nPoints) {
	    Tplus1=T_old[i];
	    Dplus = diffusion*(1.0-x[i]*x[i]);
	}

	if(i/=nPoints) {
	    Tplus1=T_old[i+1];
	    Dminus = 0.5*diffusion*((1.0-x[i+1]*x[i+1]) + (1.0-x[i]*x[i]));
	}

	T1 = T_old[i];

	if (i == 0)
	    {
	    dx1 = deltax[i];
	    dx = 0.5 * (deltax[i + 1] + deltax[i]);
	    }
	else if (i == nPoints)
	    {
	    dx1 = 0.5 * (deltax[i - 1] + deltax[i]);
	    dx = deltax[i];
	    }
	else
	    {
	    dx = 0.5 * (deltax[i] + deltax[i + 1]);
	    dx1 = 0.5 * (deltax[i - 1] + deltax[i]);
	    }

	  Fj = (Dplus*(Tplus1-T1)/dx - Dminus*(T1-Tminus1)/dx1)/(0.5*(dx1+dx));
	  T[i] = T1 + (dtLEBM/C[i])*(Q[i] +Fj);

	}

    }


