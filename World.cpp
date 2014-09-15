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
#include <algorithm>
#include <numeric>

double freeze = 273.0;
double boil = 373.0;
double sigma_SB = 5.67e-5;

World::World() :
	Body()
    {
    nPoints = 0;
     obliquity = 0.0;
     rotationPeriod = 0.0;
     winterSolstice = 0.0;
     oceanFraction = 0.0;
     landFraction = 1.0-oceanFraction;
     initialTemperature = 0.0;
     nFloat = float(nPoints);

     rho_moon = 5.0e-9; // density in kg m^-3
     rigid = 4e9; // rigidity in N m^-2 (Pa)
     Qtidal = 100.0;
     hostBody = 0;

     nPoints1 = nPoints+1;
     activateMelt = false;
     restart = false;
     tidalOn = false;

     dtLEBM = 0.0;
     diffusion = 0.0;
     logFile = 0;
     snapshotFile = 0;

    }
World::World(string namestring, string typestring, double m, double rad,
	Vector3D pos, Vector3D vel, int n, double obliq, double rot, double winter,
	double ocean, double T0,bool melt, bool start, bool tide) :
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

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9; // rigidity in N m^-2 (Pa)
    Qtidal = 100.0;
    hostBody = 0;

    nPoints1 = nPoints+1;
    activateMelt = melt;
    restart = start;
    tidalOn = tide;
    initialiseLEBM();

    }
World::World(string namestring, string typestring, double m, double rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, int n,
	double obliq, double rot, double winter, double ocean, double T0, bool melt, bool start, bool tide) :
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
    activateMelt = melt;
    restart = start;
    tidalOn = tide;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9;
    Qtidal = 100.0;
    hostBody = 0;

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

    int i;
    double dtmax = 1.0e30;
    // Set length of vectors

   lat.resize(nPoints1,0.0);
   x.resize(nPoints1,0.0);
   coslat.resize(nPoints1,0.0);
   tanlat.resize(nPoints1,0.0);
   deltax.resize(nPoints1,0.0);
   iceFraction.resize(nPoints1,0.0);
   hab.resize(nPoints1,0.0);
   infrared.resize(nPoints1,0.0);
   Q.resize(nPoints1,0.0);
   albedo.resize(nPoints1,0.0);
   insol.resize(nPoints1,0.0);
   tidal.resize(nPoints1,0.0);
   tau.resize(nPoints1,0.0);
   C.resize(nPoints1,0.0);

   // Set temperature equal to initial temperature
   T.resize(nPoints1,initialTemperature);
   T_old.resize(nPoints1, 0.0);

   // Set up latitudes etc

   double dlat = pi/nFloat;

   for (i=0; i< nPoints1; i++)
       {
       lat[i] = -pi/2.0 + i*dlat;
       x[i] = sin(lat[i]);
       coslat[i] = cos(lat[i]);
       if(coslat[i]<0.0){coslat[i] = 0.0;}
       tanlat[i] = tan(lat[i]);
       deltax[i] = coslat[i]*dlat;
       }

   // Set up diffusion constant

   diffusion = 5.394e2 * rotationPeriod*rotationPeriod;

   // Set up variables to handle ice melting
   meltTime.resize(nPoints1,0.0);
   melting.resize(nPoints1,false);

   // Set up files

   initialiseOutputVariables(restart);
   // Calculate initial parameters

#pragma omp parallel default(none) \
	shared(freeze,boil)\
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < nPoints1; i++)
	    {
	    calcIce(i);
	    calcHeatCapacity(i);
	    calcOpticalDepth(i);
	    calcAlbedo(i);
	    calcCooling(i);
	    calcNetHeating(i);

	    if(tidalOn && hostBody!=0){
		cout << "calculating heating " << endl;
		calcTidalHeating(i);}

	    calcHabitability(i,freeze,boil);
	    }
	}
   calcLEBMTimestep(dtmax);

    }

void World::setTemperature(vector<double>temp){

    int i;
    double dtmax = 1.0e30;
	for(i=0; i<nPoints; i++)
	    {
	    T[i] = temp[i];
	    }
	T[nPoints] = T[nPoints-1];

	   // Calculate initial parameters

	#pragma omp parallel default(none) \
		shared(freeze,boil)\
		private(i)
		{
	#pragma omp for schedule(runtime) ordered
		for (i = 0; i < nPoints1; i++)
		    {
		    calcIce(i);
		    calcHeatCapacity(i);
		    calcOpticalDepth(i);
		    calcAlbedo(i);
		    calcCooling(i);
		    if(tidalOn) {calcTidalHeating(i);}
		    calcNetHeating(i);
		    calcHabitability(i,freeze,boil);
		    }
		}
	   calcLEBMTimestep(dtmax);

    }


void World::updateLEBM(vector<Body*> bodies, vector<double>eclipsefrac, double &dtmax)
    {
    /*
     * Written 10/1/14 by dh4gan
     * This method is an umbrella method, using the other methods
     * written here to advance the LEBM model by dt, for convenient use by other methods
     * This assumes the timestep has already been calculated
     */
    setInsolationZero();
    int bodyCount = bodies.size();
    int i,b;

    for (b=0; b< bodyCount; b++)
	{


	if (bodies[b]->getType()=="Star" or bodies[b]->getType()=="Planet")
	    {
	    calcInsolation(bodies[b],eclipsefrac[b]);
	    }

	}

#pragma omp parallel default(none) \
	shared(freeze,boil)\
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < nPoints1; i++)
	    {
	    calcIce(i);
	    calcHeatCapacity(i);
	    calcOpticalDepth(i);
	    calcAlbedo(i);
	    calcCooling(i);

	    if(tidalOn && hostBody!=0)
		{
		calcTidalHeating(i);
		}
	    if(tidalOn && hostBody==0)
		{
		cout << "Warning: Host Body undefined, tidal heating inactive" << endl;
		}

	    calcNetHeating(i);

	    calcHabitability(i,freeze,boil);
	    }
	}

	calcLEBMTimestep(dtmax);

    integrate();

    }

void World::updateLEBM(vector<Body*> bodies, vector<double> eclipsefrac)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Overloaded method, with dtLEBM enforced from above
     */

    double dtmax = 1.0e30;

    updateLEBM(bodies,eclipsefrac, dtmax);

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

    int i;
    double cos_H;
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

    double sind = sin(declination);
    double cosd = cos(declination);
    double tand = tan(declination);

    // Insolation prefactor depends on luminosity and separation only

    double lstar = star->getLuminosity();

#pragma omp parallel default(none) \
	shared(sind,cosd,tand,fluxsolcgs)\
	shared(pi,magpos,lstar,eclipsefrac)\
	private(i,cos_H,H)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < nPoints1; i++)
	    {

	    // calculate the diurnally averaged hour angle

	    cos_H = -tanlat[i] * tand;

	    if (fabs(cos_H) > 1.0)
		cos_H = cos_H / fabs(cos_H);

	    H = acos(cos_H);

	    insol[i] = insol[i]
		    + fluxsolcgs * lstar * (1.0 - eclipsefrac)
			    / (pi * magpos * magpos)
			    * (H * x[i] * sind + coslat[i] * cosd * sin(H));
	    if(insol[i]>1.0e10){

		cout << i << star->getName() <<  "  "<<insol[i] <<"  " << fluxsolcgs * lstar * (1.0 - eclipsefrac)
				    / (pi * magpos * magpos)
				    * (H * x[i] * sind + coslat[i] * cosd * sin(H)) << "  "  <<lstar <<  "  " << (1.0 - eclipsefrac)
				    << "  " << (pi * magpos * magpos) << endl;}
	    }

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
    double C_ice = 0.0;
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

void World::calcTidalHeating(int iLatitude)
    {
    /*
     * Written 13/8/14 by dh4gan
     * Calculates the tidal heating on the World given the orbital elements
     *
     */


    // Calculate the current orbit of the World around the host
    calcOrbitFromVector(Gmau, hostBody);

    double hostMass = hostBody->getMass()*msol;

    //Calculate tidal heating in SI units

    tidal[iLatitude] = 21 * (pow(Gsi * hostMass, 2.5))
	    * pow(3 * mass*msol / (4.0 * pi), 1.666) * eccentricity * eccentricity
	    * pow(rho_moon, 0.333);

    tidal[iLatitude] = tidal[iLatitude]/ (38.0 * rigid * Qtidal * pow(semiMajorAxis * AU, 7.5));
    tidal[iLatitude] = tidal[iLatitude]*1.0e3; // convert this into cgs


    }

void World::calcNetHeating(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculates the Net Heating in the LEBM system
     *
     */

	Q[iLatitude] = insol[iLatitude]*(1.0-albedo[iLatitude]) + tidal[iLatitude]- infrared[iLatitude];
    }

void World::calcHabitability(int iLatitude, double &minT, double &maxT)
    {
    /*
     * Written 9/1/14 by dh4gan
     *
     * Calculates a habitability index based on a specified set of
     * minimum and maximum temperatures for life
     */

    if (T[iLatitude] >= minT and T[iLatitude] <= maxT)
	{
	hab[iLatitude] = 1;
	}
    else
	{
	hab[iLatitude] = 0.0;
	}

    }

void World::calcLEBMTimestep(double &dtmax)
    {
/*
 * Written 10/1/14 by dh4gan
 * Calculates the minimum timestep for the LEBM (in seconds)
 *
 */
    int i;
    vector<double> timestep(nPoints,0);
    double Dplus;

    dtLEBM = 1.0e30;

// Calculate the timestep at every latitude on the World

#pragma omp parallel default(none) \
	shared(timestep)\
	private(i,Dplus)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < nPoints1; i++)
	    {
	    if (i == nPoints)
		{
		Dplus = diffusion * (1.0 - x[i] * x[i]);
		}
	    else
		{
		Dplus = diffusion * 0.5
			* ((1.0 - x[i] * x[i]) + (1.0 - x[i + 1] * x[i + 1]));
		}

	    timestep[i] = 1.0e30;

	    if (i > 0 and i < nPoints)
		{
		timestep[i] = deltax[i] * deltax[i] * C[i] / (2.0 * Dplus);
		}
	    }
	}

// Now find the minimum timestep



    dtLEBM = *(min_element(timestep.begin(), timestep.end()));

    if(dtLEBM==1.0e30)
	{
	cout << "ERROR in LEBM Timestep calculation for World " << getName() <<  endl;
	cout << dtLEBM << "  " << diffusion << endl;
	exit(EXIT_FAILURE);
	dtLEBM = -10.0;
	}

    // If bigger than imposed maximum timestep, then set it to this value
    if(dtLEBM > dtmax){
	dtLEBM = dtmax;
    }


    }

void World::integrate()
    {
    /*
     * Written 10/1/14 by dh4gan
     * Integrates the diffusion equation to drive the system forward
     * This method includes an ice melting algorithm
     */
    int i;
    double Tminus1,Tplus1, Dplus,Dminus;
    double T1,dx,dx1,Fj, period;
    double year = 3.15e7;

    T_old = T;

	// For melting algorithm, find period by calling the Base class method
	period = Body::getPeriod();

	// Convert into seconds
	// WARNING - THIS DEPENDS ON CHOICE OF G, and mass and distance units
	// CANONICAL UNITS USED IN THIS CODE: G=1, M = Msol, d = AU ==> t = 2 pi units/year

	period = period*year/(2.0*pi);

#pragma omp parallel default(none) \
	shared(freeze,period) \
	private(i,Tminus1,Dminus,Tplus1,Dplus)\
	private(T1,dx,dx1, Fj)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < nPoints1; i++)
	    {

	    // If ice melting algorithm began in a previous timestep, modify albedo
	    // melting ice has a slightly lower albedo (Spiegel et al 2010)

	    if(melting[i])
	    	{

	    	albedo[i] = albedo[i]*0.8;
	    	}

	    // Ifs ensure we use the correct ghost cells to do the integration
	    // Do the minimum value first

	    if (i == 0)
		{
		Tminus1 = T_old[i];
		Dminus = diffusion * (1.0 - x[i] * x[i]);
		}
	    else
		{
		Tminus1 = T_old[i - 1];
		Dminus = 0.5 * diffusion
			* ((1.0 - x[i - 1] * x[i - 1]) + (1.0 - x[i] * x[i]));
		}

	    // Now the maximum value

	    if (i == nPoints)
		{
		Tplus1 = T_old[i];
		Dplus = diffusion * (1.0 - x[i] * x[i]);
		}
	    else
		{
		Tplus1 = T_old[i + 1];
		Dplus = 0.5 * diffusion
			* ((1.0 - x[i + 1] * x[i + 1]) + (1.0 - x[i] * x[i]));
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

	    Fj = (Dplus * (Tplus1 - T1) / dx - Dminus * (T1 - Tminus1) / dx1)
		    / (0.5 * (dx1 + dx));
	    T[i] = T1 + (dtLEBM / C[i]) * (Q[i] + Fj);

	    //if(T[i]>1.0e6 or T[i]<0.0){
	    //cout <<i << "  "<< T[i] << "   " << "  " << T1 <<"  " << dtLEBM << "  " << (dtLEBM / C[i]) * (Q[i] + Fj) << "  " << insol[i] << "    " << infrared[i] << "   " << Q[i] << endl;
	    //}
	    // If this temperature brings a cold world to the freezing point, then begin ice
	    // melting algorithm

	    // Need to ensure melting algorithm activated in one way only

	    if(T[i]>=freeze and T_old[i] < freeze and activateMelt)
		{
		melting[i]=true;
		//cout << "MELTING " << i << "  " << nPoints1 << endl;
		}

	    // If this latitude is currently melting and net heating is positive,
	    // then ice sheet continues to melt

	    if (melting[i])
		{
		if (Q[i] >= 0.0)
		    {
		    meltTime[i] = meltTime[i] + dtLEBM;
		    }
		else
		    {
		    meltTime[i] = 0.0;
		    }



		// If melting process continues for more than one orbit, then ice sheet has melted
		if (meltTime[i] < period)
		    {
		    T[i] = freeze-1.0e-2;
		    }

		if(meltTime[i]>=period)
		    {
		    melting[i] = false;
		    }

		}
	    // cout <<" Integrate: "<<  i <<"   "<< T[i] <<"   "<< dtLEBM<<"   " << Q[i]<<"   " << Fj <<endl;

	    }
	}

    }

void World::initialiseOutputVariables(bool restart)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Sets up the log file for this World
     *
     */

    string logFileName = getName()+".log";

    if(restart)
	{
	logFile = fopen(logFileName.c_str(), "a");
	}
    else
	{
	logFile = fopen(logFileName.c_str(), "w");
	}
    }

void World::calcLEBMMeans(double &minT, double &maxT, double &meanT, double &meanQ, double &meanA, double &meanIR, double &meanS,
	double &meanhab, double &meanTidal)
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
    meanTidal = 0.0;

    for (int i=0; i< nPoints1; i++)
	{
	meanT = meanT + T[i]*0.5*deltax[i];
	meanQ = meanQ + Q[i]*0.5*deltax[i];
	meanA = meanA + albedo[i]*0.5*deltax[i];
	meanIR = meanIR + infrared[i]*0.5*deltax[i];
	meanS = meanS + insol[i]*0.5*deltax[i];
	meanhab = meanhab + hab[i]*0.5*deltax[i];
	meanTidal = meanTidal + tidal[i]*0.5*deltax[i];

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

    double minT, maxT, meanT, meanQ, meanA, meanIR, meanS, meanhab, meanTidal;


    // Firstly, write line to log file
    calcLEBMMeans(minT, maxT, meanT, meanQ, meanA, meanIR,meanS, meanhab, meanTidal);

    // Also include orbital data here
    fprintf(logFile, "%+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+4.E \n",
	    tSnap, minT, maxT, meanT,meanQ,meanA,meanIR,meanS,meanhab, meanTidal,
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
	fprintf(snapshotFile, "%+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E %+.4E \n",
		x[i], lat[i], T[i], C[i], Q[i], infrared[i],
		albedo[i], insol[i], tau[i],iceFraction[i], hab[i], tidal[i]);
	}
    fclose(snapshotFile);

    }

int World::findRestartTemperature()
    {

    /*
     * written by dh4gan 13/5/14
     * Reads last temperature snapshot file for restart purposes
     *
     */

    // Read log file to determine how many snapshots there have been

    string logFileName = getName()+".log";
    string line;
    int nSnap = 0;
    double blank;
    double dtmax = 1.0e30;
    ifstream lfile(logFileName.c_str());

    if(lfile)
   	{

   	while(getline(lfile,line))
   	    {
   	    nSnap++;
   	    }

   	cout << "Number of Snapshots is " << nSnap << endl;

   	}
       else
   	{
   	cout << "Error: File "<< logFileName << "not found" << endl;
   	return -1;
   	}

    lfile.close();

    ostringstream convert;
    convert << nSnap;

    string numString = convert.str();
    string worldFile = getName()+"."+numString;

    ifstream wfile(worldFile.c_str());

    // Read header (number of grid points, time)
    getline(wfile,line);

    int i=0;

    while (getline(wfile, line))
	{
	istringstream iss(line);
	iss >> blank;
	iss >> blank;
	if (i < nPoints1)
	    {
	    iss >> T[i];
	    }
	i++;
	}

#pragma omp parallel default(none) \
	shared(freeze,boil)\
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < nPoints1; i++)
	    {
	    calcIce(i);
	    calcHeatCapacity(i);
	    calcOpticalDepth(i);
	    calcAlbedo(i);
	    calcCooling(i);
	    calcNetHeating(i);
	    calcHabitability(i,freeze,boil);
	    }
	}
   calcLEBMTimestep(dtmax);

   return nSnap;
    }



