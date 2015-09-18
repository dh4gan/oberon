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
    type = "World";
    nPoints = 0;
     obliquity = 0.0;
     precession = 0.0;
     ellipticity = 0.00328005;
     rotationPeriod = 0.0;
     oceanFraction = 0.0;
     landFraction = 1.0-oceanFraction;
     initialTemperature = 0.0;
     nFloat = float(nPoints);

     luminosity = 0.0;

     rho_moon = 5.0e-9; // density in kg m^-3
     rigid = 4e9; // rigidity in N m^-2 (Pa)
     Qtidal = 100.0;


     nPoints1 = nPoints+1;
     activateMelt = false;
     restart = false;

     tidalHeatingOn = false;
     obliquityEvolutionOn = false;

     CScycleOn = false;

     dtLEBM = 0.0;
     diffusion = 0.0;
     logFile = 0;
     latFile = 0;
     snapshotFile = 0;

    }
World::World(string namestring, double m, double rad,Vector3D pos, Vector3D vel, int n, double obliq, double rot, double prec,
	double ocean, double T0,bool melt, bool start, bool tide, bool obevol, bool CScycle) :

	Body(namestring, m, rad, pos, vel)
    {
    type = "World";
    nPoints = n;
    obliquity = obliq;
    ellipticity = 0.00328005;
    precession = prec;
    rotationPeriod = rot;
    oceanFraction = ocean;
    landFraction = 1.0-oceanFraction;
    initialTemperature = T0;
    nFloat = float(nPoints);

    luminosity = 0.0;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9; // rigidity in N m^-2 (Pa)
    Qtidal = 100.0;


    nPoints1 = nPoints+1;
    activateMelt = melt;
    restart = start;

    tidalHeatingOn = tide;
    obliquityEvolutionOn = obevol;

    CScycleOn = CScycle;
    initialiseLEBM();

    }
World::World(string namestring, double m, double rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, int n,

	double obliq, double rot, double prec, double ocean, double T0, bool melt, bool start, bool tide, bool obevol, bool CScycle) :
	Body(namestring, m, rad, semimaj, ecc, inc, longascend,
		argper, meananom, G, totalMass)
    {
    type = "World";
    nPoints = n;
    obliquity = obliq;
    ellipticity = 0.00328005;
    precession = prec;
    rotationPeriod = rot;
    oceanFraction = ocean;
    landFraction = 1.0-oceanFraction;
    initialTemperature = T0;
    nFloat = float(nPoints);
    nPoints1 = nPoints+1;
    activateMelt = melt;
    restart = start;

    tidalHeatingOn = tide;
    obliquityEvolutionOn = obevol;

    luminosity = 0.0;
    CScycleOn = CScycle;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9;
    Qtidal = 100.0;

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
   //CScycle.resize(nPoints1,0.0); //Not sure if necessary... Giblin 10/7/15
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

   calcLuminosity();

   // Set up files

   initialiseOutputVariables(restart);
   // Calculate initial parameters

#pragma omp parallel default(none) \
	shared(freeze,boil,cout)\
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

	    if(tidalHeatingOn && hostBody!=0){
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
		    if(tidalHeatingOn) {calcTidalHeating(i);}
		    calcNetHeating(i);
		    calcHabitability(i,freeze,boil);
		    }
		}
		calcLuminosity();
	   calcLEBMTimestep(dtmax);

    }


void World::calcObliquity(vector<Body*> bodies, double G, double totmass)
    {

    /*
     * Written 01/04/16 by dh4gan
     * Calculates the obliquity and precession evolution
     * of the World using Laskar et al (1986a,b) functions
     * (see also Armstrong et al 2014)
     *
     */


    // Store orbital data from previous timestep
    double oldInclination = inclination;
    double oldAscending = longitudeAscendingNode;

    // Calculate new orbital parameters

    if(hostBody){
	calcOrbitFromVector(G,hostBody);
    }
    else
	{
    calcOrbitFromVector(G,totmass);
	}

    // Calculate rate of change of inclination and ascending node

    double dinc = inclination-oldInclination;

    // Be careful with cases where inclination crosses 0->2pi line
    if(dinc < -pi) dinc = dinc+twopi;
    if(dinc > pi) dinc = dinc-twopi;


    dinc = dinc*year/dtLEBM;

    double domega = longitudeAscendingNode - oldAscending;

    // Again be careful with 0-> 2pi line

    if(domega < -pi) domega = domega+twopi;
    if(domega > pi) domega = domega-twopi;

    domega = domega*year/dtLEBM;

   if(fabs(dinc) > 0.5*year/dtLEBM)
{
    cout << "inc: " << inclination << "   " << oldInclination << "   " << dinc*dtLEBM/year << endl;
}
if(fabs(domega)>0.5*year/dtLEBM)
{ 
   cout << "ome: " << longitudeAscendingNode << "   " << oldAscending << "   " << domega*dtLEBM/year << endl;
}
    // Calculate p and q parameters

    double p = sin(0.5*inclination)*sin(longitudeAscendingNode);
    double q = sin(0.5*inclination)*cos(longitudeAscendingNode);

    double pdot = 0.5*cos(0.5*inclination)*sin(longitudeAscendingNode)*dinc + q*domega;
    double qdot = 0.5*cos(0.5*inclination)*cos(longitudeAscendingNode)*dinc -p*domega;

    // A,B and C functions

    double Cfunc = q*pdot - p*qdot;
    double prefactor =2.0/(sqrt(1.0-p*p-q*q));
    double Afunc = prefactor*(qdot -p*Cfunc);
    double Bfunc = prefactor*(pdot -q*Cfunc);

    // Direct torques from the host body
    // If no host, then calculate direct torques and relativistic precessions for all stars in the system
    // ASSUMING SINGLE ORBITAL SOLUTION (e.g. around centre of mass) APPROPRIATE FOR EACH STAR

    double minuse2 = 1.0 - eccentricity*eccentricity;
    double Rtorque = 0.0;
    double pGR = 0.0;

    if(hostBody)
	{
	double bodyMass = hostBody->getMass();
	double k_Kep = Gmau_day*bodyMass/(twopi*twopi);
	Rtorque= 3.0*k_Kep*rotationPeriod/(twopi*semiMajorAxis*semiMajorAxis*semiMajorAxis);
	Rtorque = Rtorque*ellipticity*cos(obliquity)*(0.5*pow(minuse2,-1.5)-0.522e-6);

	pGR = pow(k_Kep,1.5)*pow(1.0+mass/bodyMass,0.5)/(2.0*minuse2*c_mau*c_mau*pow(semiMajorAxis,2.5));


	}
    else
	{
	// Find all stars in the System, and compute torques and relativistic precession
	// ASSUMING SINGLE ORBITAL SOLUTION APPROPRIATE FOR EACH STAR

	for (int ibody=0; ibody < int(bodies.size()); ibody++)
	    {
	    if(bodies[ibody]->getType()=="Star")
		{
		double bodyMass = bodies[ibody]->getMass();
		double k_Kep = G*bodyMass/(twopi*twopi);
		double Rtorquestar= 3.0*k_Kep*rotationPeriod/(twopi*semiMajorAxis*semiMajorAxis*semiMajorAxis);
		Rtorquestar = Rtorquestar*ellipticity*cos(obliquity)*(0.5*pow(minuse2,-1.5)-0.522e-6);
		Rtorque = Rtorque + Rtorquestar*0.002737;  // Correction factor comes from using years, not days

		pGR = pGR + pow(k_Kep,1.5)*pow(1.0+mass/bodyMass,0.5)/(2.0*minuse2*c_mau*c_mau*pow(semiMajorAxis,2.5));

		}
	    }

	}

    // Now update obliquity and precession parameters
    // (Note: cot(x) replaced by cos(x)/sin(x)

    double cotObliq = 0.0;
    if(obliquity>0.0){
	cotObliq = cos(obliquity)/sin(obliquity);
    }


    double precdot = Rtorque - cotObliq*(Afunc*sin(precession)+Bfunc*cos(precession)) -2.0*Cfunc - pGR;
    double obliqdot = -Bfunc*sin(precession) +Afunc*cos(precession);

    precession = precession + precdot*dtLEBM/year;
    obliquity = obliquity + obliqdot*dtLEBM/year;


 //   cout  << "p, q: " << p << "  " << q << "   " << pdot << "  " << qdot << endl;
//    cout << "OBLIQUITY: " <<obliquity << "   " <<  obliqdot<<"  " << Bfunc <<  "  " << Afunc << "  " << precession << endl;
    precession = fmod(precession,twopi);

    if(precession < 0.0){precession = twopi+precession;}

    obliquity = fmod(obliquity,twopi);
    if(obliquity < 0.0){obliquity = twopi+precession;}
    //cout << "PRECESSION " << precession << endl;


    }

void World::calcLuminosity()
    {
/*
 * Written 11/8/14 by dh4gan
 * Calculates the luminosity of the Planet, given its temperature and reflected starlight
 *
 */
    double pi = 3.141592653;
    double sigma_SB =5.67e-8;
    double AU = 1.496e11;
    double lsol = 3.8626e26;

    double minT, maxT,meanT,meanQ,meanA,meanIR,meanS,meanhab,meanTidal;

    calcLEBMMeans(minT, maxT, meanT, meanQ, meanA, meanIR,meanS, meanhab, meanTidal);

    luminosity = 4.0*pi*getRadius()*getRadius()*AU*AU*sigma_SB*meanT*meanT*meanT*meanT/lsol;

    }

void World::updateLEBM(vector<Body*> bodies, double &G, double &totmass, vector<double>eclipsefrac, double &dtmax, bool &planetaryIllumination)
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

    // Compute current obliquity and precession vector
    if(obliquityEvolutionOn)
	{
    calcObliquity(bodies,G,totmass);
	}

    for (b=0; b< bodyCount; b++)
	{

    	// Skip if trying to calculate self-insolation
    if(bodies[b]->getName() == getName() ) continue;

	if (bodies[b]->getType()=="Star")
	    {
	    calcInsolation(bodies[b],eclipsefrac[b]);
	    }

	if((bodies[b]->getType()=="Planet" or bodies[b]->getType()=="World") and (planetaryIllumination==true))
	{
		    calcInsolation(bodies[b],eclipsefrac[b]);
	}

	}

#pragma omp parallel default(none) \
	shared(freeze,boil,cout)\
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

	    if(tidalHeatingOn && hostBody!=0)
		{
		calcTidalHeating(i);
		}
	    if(tidalHeatingOn && hostBody==0)
		{
		cout << "Warning: Host Body undefined, tidal heating inactive" << endl;
		}

	    calcNetHeating(i);

	    calcHabitability(i,freeze,boil);
	    }
	}

	calcLuminosity();
	calcLEBMTimestep(dtmax);

    integrate();

    }

void World::updateLEBM(vector<Body*> bodies, double &G, double &totmass, vector<double> eclipsefrac,bool &planetaryIllumination)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Overloaded method, with dtLEBM enforced from above
     */

    double dtmax = 1.0e30;

    updateLEBM(bodies,G,totmass,eclipsefrac, dtmax,planetaryIllumination);

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

    Vector3D decVector = orbitalAngularMomentum.unitVector();

    // Rotate this vector if world has non-zero obliquity
    if (obliquity > 0.0) {
	decVector.rotateY(-obliquity);
	}


    // rotate decVector to account for precession
    // This is a rotation of the spin-axis around the orbital angular momentum axis (not Z)

    decVector.rotateAboutAxis(orbitalAngularMomentum, precession);

    // Obtain declination angle

    double rdotn = unitpos.dotProduct(decVector);
    //Vector3D rcrossn = unitpos.crossProduct(decVector);
    //double rcrossnMag = rcrossn.dotProduct(rcrossn.unitVector());


    double declination = safeAcos(rdotn);


    declination = piby2-declination;
     //cout << "ROTATIONS: " << endl;
    //orbitalAngularMomentum.printVector();
    //decVector.printVector();

    // Allowed range of declinations: -pi, pi
    // cos is an even function!


    if(declination > obliquity) declination = pi-declination;

    //if (decVector.elements[1] <0.0) declination = -1*declination;

    double sind = sin(declination);
    double cosd = cos(declination);
    double tand = tan(declination);


    //cout <<"Declination: " << trueAnomaly << "   " <<  obliquity << "   " << rdotn << "   " << rcrossnMag << "   "<< declination << "   " << decCross << endl;

    // Insolation prefactor depends on luminosity and separation only

    double lstar = star->getLuminosity();

#pragma omp parallel default(none) \
	shared(sind,cosd,tand,fluxsolcgs,star,cout)\
	shared(magpos,lstar)\
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

		if(insol[i]<0.0) {insol[i]=0.0;}
	    if(insol[i]>1.0e10 or insol[i]!=insol[i] or insol[i] < 0.0){

	      cout << "ERROR: Negative insolation calculated " << endl;
	      cout << i << "  " << star->getName () << "  " << insol[i] << "  "
		  << precession << "  " << obliquity << "   " << "   " << cosd << "  "
		  << sind << "  " << tand << endl;


	    }
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









void World::calcC02Pressure(int iLatitude)
    {
    /*
     * Written 10/7/15 by BenjaminGiblin
     * Calculates C02 pressure using regime defined
     * on pg. 5 of Spiegel et al 2010
     */

	if (T[iLatitude] >= 290.)
	{
	CO2pressure[iLatitude] = 3.3e-4; //bars
	}
    else if (T[iLatitude] < 290. and T[iLatitude] > 250.)
	{
	CO2pressure[iLatitude] = pow(10., -2-(T[iLatitude]-250.)/27.); //bars
	}

    else if (T[iLatitude] <= 250.)
	{
	CO2pressure[iLatitude] = 0.01; //bars
	}
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
	cout << "Timesteps: " << endl;
	for(i=0; i<nPoints1; i++)
	{
		cout << timestep[i] << "   " << T[i] << "   " << C[i] << endl;
	}

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
	shared(freeze,period,cout) \
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

	    if(T[i]>1.0e6 or T[i]<0.0){
	    cout << "Temperature out of sensible range " << endl;
	    cout <<i << "  "<< T[i] << "   " << "  " << T1 <<"  " << dtLEBM << "  " << (dtLEBM / C[i]) * (Q[i] + Fj) << "  " << insol[i] << "    " << infrared[i] << "   " << Q[i] << endl;
	    }

	    // This statement checks for NaNs
	    if(T[i]!=T[i]){
	    		cout << "NaN detected in T " << endl;
	    	    cout <<i << "  "<< T[i] << "   " << "  " << T1 <<"  " << dtLEBM << "  " << C[i] << "  " << Q[i] << "  " <<Fj  << "  " << insol[i] << "    " << infrared[i] << "   " << Q[i] << endl;
	    	    exit(EXIT_FAILURE);
	    	    }
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
    string latFileName = getName()+".lat";

    if(restart)
	{
	logFile = fopen(logFileName.c_str(), "a");
	latFile = fopen(latFileName.c_str(), "a");
	}
    else
	{
	logFile = fopen(logFileName.c_str(), "w");
	latFile = fopen(latFileName.c_str(), "w");

	// Write header for latitudinal file

	for (int i=0; i<nPoints; i++){
	fprintf(latFile, "%+.6E  ", lat[i]);
	}

	fprintf(latFile, "\n");
	fflush(latFile);

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

void World::outputLEBMData(int &snapshotNumber, double &tSnap, bool fullOutput)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Outputs a snapshot of the LEBM to file
     * and adds a line to the log file
     *
     */

    double minT, maxT, meanT, meanQ, meanA, meanIR, meanS, meanhab, meanTidal;
    string formatString;

    // Firstly, write line to log file
    calcLEBMMeans(minT, maxT, meanT, meanQ, meanA, meanIR,meanS, meanhab, meanTidal);

    // Also include orbital data here
    formatString = "+%.6E  ";

    for (int icol=0;icol < 17; icol++)
	{
	formatString = formatString + "%.6E  ";
	}
    formatString = formatString + "\n";

    fprintf(logFile, formatString.c_str(),
	    tSnap, minT, maxT, meanT,meanQ,meanA,meanIR,meanS,meanhab, meanTidal,
	    semiMajorAxis, eccentricity, inclination,
	    argumentPeriapsis, longitudeAscendingNode, meanAnomaly,
	    obliquity*radToDeg,precession*radToDeg);
    fflush(logFile);

    // Write latitudinal temperature data

    fprintf(latFile, "%+.6E  ", tSnap);

    for (int i=0; i<nPoints; i++){
	fprintf(latFile, "%+.6E  ", T[i]);
    }

    fprintf(latFile, "\n");
    fflush(latFile);

    if(fullOutput){

    // If going for full LEBM output:
    // Now write snapshot of LEBM

    ostringstream convert;
    convert << snapshotNumber;

    string numString = convert.str();
    string snapshotFileName = getName()+"."+numString;

    snapshotFile = fopen(snapshotFileName.c_str(), "w");

    fprintf(snapshotFile, "%i %+.6E \n", nPoints, tSnap);


    for (int i=0; i<nPoints; i++)
	{
	fprintf(snapshotFile, "%+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E \n",
		x[i], lat[i], T[i], C[i], Q[i], infrared[i],
		albedo[i], insol[i], tau[i],iceFraction[i], hab[i], tidal[i]);
	}
    fclose(snapshotFile);
    }

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



