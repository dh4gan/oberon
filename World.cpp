/*
 * World.cpp
 *
 *  Created on: Jan 9 2014
 *      Author: dh4gan
 */

#include "World.h"
#include "Constants.h"
#include <stdlib.h>
#include <fstream>
#include <sstream>


double freeze = 273.0;
double boil = 373.0;
double sigma_SB = 5.67e-5;


World::World() :
	Body()
    {
    }
World::World(string namestring, string typestring, double m, double rad,
	Vector3D pos, Vector3D vel, int n, double obliq, double rot, double winter,
	double ocean, double T0) :
	Body(namestring, typestring, m, rad, pos, vel)
    {
    nPoints = n;
    obliquity = obliq;
    rotationPeriod = rot;
    winterSolstice = winter;
    oceanFraction = ocean;
    landFraction = 1.0-oceanFraction;
    initialTemperature = T0;
    nFloat = float(nPoints);

    nPoints1 = nPoints+1;

    initialiseLEBM();

    }
World::World(string namestring, string typestring, double m, double rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, int n,
	double obliq, double rot, double winter, double ocean, double T0) :
	Body(namestring, typestring, m, rad, semimaj, ecc, inc, longascend,
		argper, meananom, G, totalMass)
    {
    nPoints = n;
    obliquity = obliq;
    rotationPeriod = rot;
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

   // Set up files

   initialiseOutputVariables();
   // Calculate initial parameters

   for (int i=0; i<nPoints1; i++)
       {
	calcIce(i);
	calcHeatCapacity(i);
	calcOpticalDepth(i);
	calcAlbedo(i);
	calcCooling(i);
	calcNetHeating(i);
       }

   calcHabitability(freeze,boil);
   calcLEBMTimestep();

    }


void World::updateLEBM(vector<Body*> bodies, vector<double>eclipsefrac)
    {
    /*
     * Written 10/1/14 by dh4gan
     * This method is an umbrella method, using the other methods
     * written here to advance the LEBM model by dt, for convenient use by other methods
     * This assumes the timestep has already been calculated
     */
    setInsolationZero();
    int bodyCount = bodies.size();

    for (int b=0; b< bodyCount; b++)
	{

	if (bodies[b]->getType()=="Star")
	    {
	    calcInsolation(bodies[b],eclipsefrac[b]);
	    }
	}

    for (int i = 0; i < nPoints1; i++)
	{
	calcIce(i);
	calcHeatCapacity(i);
	calcOpticalDepth(i);
	calcAlbedo(i);
	calcCooling(i);
	calcNetHeating(i);
	}

    calcHabitability(freeze,boil);
    integrate();

    }

void World::updateLEBM(vector<Body*> bodies, vector<double> eclipsefrac, double dt)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Overloaded method, with dtLEBM enforced from above
     */

    dtLEBM = dt;
    updateLEBM(bodies,eclipsefrac);

    }


void World::calcInsolation(Body* star, double &eclipsefrac)
    {

    /*
     * Written by dh4gan, 10/1/14
     * Calculates the flux received at each latitude from input star
     * Given that eclipsefrac of the stellar light is blocked by other objects
     * This is the only calculation that has its own for loop in the updateLEBM algorithm
     * Mainly because there is a good deal of pre-calculation required to determine declinations
     * etc
     */

    vector<double> cos_H(nPoints1,0.0);
    double H;

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
	decVector.rotateX(-obliquity);
	}

    // TODO - Check implementation Rotate in y-axis if world's winter solstice longitude non-zero

    //if(winterSolstice!=0.0)
//	{
//	decVector.rotateY(winterSolstice);
//	}

    // Obtain declination angle

    double rdotn = unitpos.dotProduct(decVector);
    double declination = safeAcos(rdotn);

    // Check: is declination greater than pi?

    if (decVector.elements[1] <0.0) declination = -1*declination;
   // cout << decVector.elements[0] << "   " << decVector.elements[1] <<"   " << decVector.elements[2] << endl;


    double sind = sin(declination);
    double cosd = cos(declination);
    double tand = tan(declination);

    // Insolation prefactor depends on luminosity and separation only

    double lstar = star->getLuminosity();

    for (int i=0; i<nPoints1; i++)
	{

	// calculate the diurnally averaged hour angle

	cos_H[i] = -tan(lat[i])*tand;

	if(fabs(cos_H[i]) >1.0) cos_H[i] = cos_H[i]/fabs(cos_H[i]);

	H = acos(cos_H[i]);

	insol[i] = insol[i]+fluxsolcgs*lstar*(1.0-eclipsefrac)/(pi*magpos*magpos)*(H*x[i]*sind + cos(lat[i])*cosd*sin(H));
	//cout <<i << "   " << insol[i] << "  " <<cos_H[i] <<  "   " <<lat[i] << "  " << tan(lat[i]) << "  " <<  tand << endl;
	}


    }


void World::calcAlbedo(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Simple calculation of the albedo as a function of latitude according to Temperature
     *
     */

    albedo[iLatitude] = 0.525 - 0.245 * tanh((T[iLatitude] - freeze + 5.0) / 5.0);

    }

void World::calcHeatCapacity(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculation of Heat Capacity of the atmosphere as a function of the local temperature
     *
     */
    double C_ice;
    double C_land = 5.25e9;
    double C_ocean = 40.0 * C_land;

    if (T[iLatitude] >= freeze)
	{
	C_ice = 0.0;
	}
    else if (T[iLatitude] < freeze and T[iLatitude] > freeze - 10.0)
	{
	C_ice = 9.2 * C_land;
	}

    else if (T[iLatitude] <= freeze - 10)
	{
	C_ice = 2.0 * C_land;
	}

    C[iLatitude] = landFraction * C_land
	    + oceanFraction
		    * (iceFraction[iLatitude] * C_ice
			    + (1.0 - iceFraction[iLatitude]) * C_ocean);

    }

void World:: calcIce(int iLatitude)
    {

    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Ice Fraction given the temperature T
     *
     */

	iceFraction[iLatitude] = 1.0 - exp(-(freeze-T[iLatitude])/10.0);
	if (iceFraction[iLatitude] <0.0) iceFraction[iLatitude]=0.0;

    }

void World::calcOpticalDepth(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the optical depth of the atmosphere
     */

	tau[iLatitude] = 0.79*pow(T[iLatitude]/freeze,3);
    }

void World::calcCooling(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Infrared Cooling as a function of optical depth and temperature
     *
     */
	infrared[iLatitude] = sigma_SB*pow(T[iLatitude],4)/(1.0+0.75*tau[iLatitude]);
    }

void World::calcNetHeating(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Net Heating in the LEBM system
     *
     */

	Q[iLatitude] = insol[iLatitude]*(1.0-albedo[iLatitude]) - infrared[iLatitude];
    }

void World::calcHabitability(double &minT, double &maxT)
    {
    /*
     * Written 9/1/14 by dh4gan
     *
     * Calculates a habitability index based on a specified set of
     * minimum and maximum temperatures for life
     */

    for (int i=0; i< nPoints1; i++)
	{
	if(T[i]>=minT and T[i]<=maxT)
	    {hab[i]=1;}
	else
	    {hab[i]=0.0;}
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
	cout << "ERROR in LEBM Timestep calculation for World " << getName() <<  endl;
	cout << dtLEBM << "  " << diffusion << endl;
	exit(EXIT_FAILURE);
	dtLEBM = -10.0;
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

	// Ifs ensure we use the correct ghost cells to do the integration
	// Do the minimum value firsty

	if(i==0){
	    Tminus1 = T_old[i];
	    Dminus = diffusion*(1.0-x[i]*x[i]);
	}
	else{
	    Tminus1 = T_old[i-1];
	    Dminus = 0.5*diffusion*((1.0-x[i-1]*x[i-1]) + (1.0-x[i]*x[i]));
	}

	// Now the maximum value

	if(i==nPoints) {
	    Tplus1=T_old[i];
	    Dplus = diffusion*(1.0-x[i]*x[i]);
	}
	else {
	    Tplus1=T_old[i+1];
	    Dplus = 0.5*diffusion*((1.0-x[i+1]*x[i+1]) + (1.0-x[i]*x[i]));
	}

	// Now do the integration
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

	 // cout <<" Integrate: "<<  i <<"   "<< T[i] <<"   "<< dtLEBM<<"   " << Q[i]<<"   " << Fj <<endl;

	}

    }

void World::initialiseOutputVariables()
    {
    /*
     * Written 10/1/14 by dh4gan
     * Sets up the log file for this World
     *
     */

    string logFileName = getName()+".log";
    logFile = fopen(logFileName.c_str(), "w");
    }

void World::calcLEBMMeans(double &minT, double &maxT, double &meanT, double &meanQ, double &meanA, double &meanIR, double &meanS,
	double &meanhab)
    {
    /*
     * Written 13/1/14 by dh4gan
     * Calculates Latitudinally averaged means
     *
     */

    minT = 1.0e30;
    maxT = -1.0e30;
    meanT = 0.0;
    meanQ = 0.0;
    meanA = 0.0;
    meanIR = 0.0;
    meanS = 0.0;
    meanhab = 0.0;

    for (int i=0; i< nPoints1; i++)
	{
	meanT = meanT + T[i]*0.5*deltax[i];
	meanQ = meanQ + Q[i]*0.5*deltax[i];
	meanA = meanA + albedo[i]*0.5*deltax[i];
	meanIR = meanIR + infrared[i]*0.5*deltax[i];
	meanS = meanS + insol[i]*0.5*deltax[i];
	meanhab = meanhab + hab[i]*0.5*deltax[i];

	if(T[i] < minT) minT = T[i];
	if(T[i] > maxT) maxT = T[i];
	}


    }

void World::outputLEBMData(int &snapshotNumber, double &tSnap)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Outputs a snapshot of the LEBM to file
     * and adds a line to the log file
     *
     */

    double minT, maxT, meanT, meanQ, meanA, meanIR, meanS, meanhab;


    // Firstly, write line to log file
    calcLEBMMeans(minT, maxT, meanT, meanQ, meanA, meanIR,meanS, meanhab);

    // Also include orbital data here
    fprintf(logFile, "%+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+4.E \n",
	    tSnap, minT, maxT, meanT,meanQ,meanA,meanIR,meanS,meanhab,
	    getSemiMajorAxis(), getEccentricity(), getInclination(),
	    getArgumentPeriapsis(), getLongitudeAscendingNode(), getMeanAnomaly());
    fflush(logFile);


    // Now write snapshot of LEBM

    ostringstream convert;
    convert << snapshotNumber;

    string numString = convert.str();
    string snapshotFileName = getName()+"."+numString;
    snapshotFile = fopen(snapshotFileName.c_str(), "w");


    fprintf(snapshotFile, "%i %+.4E \n", nPoints, tSnap);

    for (int i=0; i<nPoints; i++)
	{
	fprintf(snapshotFile, "%+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E \n",
		x[i], lat[i], T[i], C[i], Q[i], infrared[i],
		albedo[i], insol[i], tau[i],iceFraction[i], hab[i]);
	}
    fclose(snapshotFile);

    }

