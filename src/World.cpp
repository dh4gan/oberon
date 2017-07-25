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
double sigma_SB = 5.67e-5; //erg s^-1 cm^-2 K^-4

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

     CSCycleOn = false;

     dtLEBM = 0.0;
     diffusion0 = 0.0;
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

    CSCycleOn = CScycle;

    outgassingRate = outgassingRateEarth;
    betaCO2 = 0.5;
    W0 = 0.0;
    gammaCO2 = 1.0;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9;
    Qtidal = 100.0;

    initialiseLEBM();

    }

World::World(string namestring, double m, double rad,Vector3D pos, Vector3D vel, int n, double obliq, double rot, double prec,
	double ocean, double T0,bool melt, bool start, bool tide, bool obevol, bool CScycle,double outgas, double beta, double seaweather, double gamma) :

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

    CSCycleOn = CScycle;

    outgassingRate = outgas;
    betaCO2 = beta;
    W0 = seaweather/year;
    gammaCO2 = gamma;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9;
    Qtidal = 100.0;

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
    CSCycleOn = CScycle;

    outgassingRate = outgassingRateEarth;
    betaCO2 = 0.5;
    W0 = 0.0;
    gammaCO2 = 1.0;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9;
    Qtidal = 100.0;

    initialiseLEBM();

    }

World::World(string namestring, double m, double rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, int n,
	double obliq, double rot, double prec, double ocean, double T0, bool melt, bool start, bool tide, bool obevol,
	bool CScycle, double outgas, double beta, double seaweather, double gamma):
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
    CSCycleOn = CScycle;

    outgassingRate = outgas;
    betaCO2 = beta;
    W0 = seaweather/year;
    gammaCO2 = gamma;

    rho_moon = 5.0e-9; // density in kg m^-3
    rigid = 4e9;
    Qtidal = 100.0;

    initialiseLEBM();

    }


World::~World()
    {
    }


void World:: setInsolationZero()
    {
    insol.assign(nPoints1,0.0);
    absorbedInsol.assign(nPoints1,0.0);
    }


void World:: resetMeanAlbedo()
    {
    /*
     * Written 13/7/2017 by dh4gan
     * Resets the arrays storing albedo averaged over all stars
     *
     */
	meanAlbedo.assign(nPoints1,0.0);
	albedoCount = 0;
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

    printf("Initialising EBM \n");
   // Set up diffusion constant
   diffusion0 = 5.394e2 * rotationPeriod*rotationPeriod;

   lat.resize(nPoints1,0.0);
   x.resize(nPoints1,0.0);
   meanZenith.resize(nPoints1,0.0);
   meanAlbedo.resize(nPoints1,0.0);
   coslat.resize(nPoints1,0.0);
   tanlat.resize(nPoints1,0.0);
   deltax.resize(nPoints1,0.0);
   iceFraction.resize(nPoints1,0.0);
   hab.resize(nPoints1,0.0);
   infrared.resize(nPoints1,0.0);
   Q.resize(nPoints1,0.0);
   albedo.resize(nPoints1,0.0);
   surfaceAlbedo.resize(nPoints1,0.0);
   insol.resize(nPoints1,0.0);
   absorbedInsol.resize(nPoints1,0.0);
   tidal.resize(nPoints1,0.0);
   CO2pressure.resize(nPoints1,CO2Earth);
   diffusion.resize(nPoints1,diffusion0);
   tau.resize(nPoints1,0.0);
   C.resize(nPoints1,0.0);

   landWeathering.resize(nPoints1,0.0);
   oceanWeathering.resize(nPoints1,0.0);
   CO2dot.resize(nPoints1,0.0);


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



   //erg cm^-2 s^-1 K^-1
   //rotationPeriod is in units of Earth's rate



   // Set up variables to handle ice melting
   meltTime.resize(nPoints1,0.0);
   melting.resize(nPoints1,false);

   resetMeanAlbedo();
   calcLuminosity();

   printf("Initialising EBM \n");
   // Set up files

   initialiseOutputVariables(restart);
   // Calculate initial parameters

   printf("Initialising EBM3 \n");
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
	    calcCooling(i);
	    calcNetHeating(i);

	    if(tidalHeatingOn && hostBody!=0){
		cout << "calculating heating " << endl;
		calcTidalHeating(i);}

	    if(CSCycleOn) {calcCO2Rates(i);}

	    calcHabitability(i,freeze,boil);
	    }
	}
	printf("Initialising EBM \n");
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
		    //calcAlbedo(i);
		    calcCooling(i);

		    if(tidalHeatingOn) {calcTidalHeating(i);}
		    if(CSCycleOn) {calcCO2Rates(i);}
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

    double minT, maxT,meanT,meanQ,meanA,meanIR,meanS,meanhab,meanTidal, meanCO2p, meanDiffusion;

    calcLEBMMeans(minT, maxT, meanT, meanQ, meanA, meanIR,meanS, meanhab, meanTidal, meanCO2p, meanDiffusion);

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
    resetMeanAlbedo();
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
	    //calcAlbedo(i);
	    calcCooling(i);

	    if(tidalHeatingOn && hostBody!=0)
		{
		calcTidalHeating(i);
		}
	    if(tidalHeatingOn && hostBody==0)
		{
		cout << "Warning: Host Body undefined, tidal heating inactive" << endl;
		}
	    if(CSCycleOn) {calcCO2Rates(i);}

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

	    meanZenith[i] = H*x[i] * sind + coslat[i] * cosd * sin(H);

	    insol[i] = insol[i]
		    + fluxsolcgs * lstar * (1.0 - eclipsefrac)
			    / (pi * magpos * magpos)
			    * (meanZenith[i]);

	    if(H>0.0)
		{
	    meanZenith[i] = meanZenith[i]/H;
		}
	    else
		{
		meanZenith[i] = 0.0;
		}

		if(insol[i]<0.0) {insol[i]=0.0;}
	    if(insol[i]>1.0e10 or insol[i]!=insol[i] or insol[i] < 0.0){

	      cout << "ERROR: Negative insolation calculated " << endl;
	      cout << i << "  " << star->getName () << "  " << insol[i] << "  "
		  << precession << "  " << obliquity << "   " << "   " << cosd << "  "
		  << sind << "  " << tand << "  " << cos_H << "   " << meanZenith[i] << endl;
	    }


	      // Note: albedo variable (and meanZenith) will get overwritten for multiple star runs!
	      // Need to calculate and use it here

	      calcAlbedo(star,i,meanZenith[i]);

	      // Also compute absorbed insolation S(1-A) here (makes net heating calculations easier)
	      absorbedInsol[i] = absorbedInsol[i]+ insol[i]*(1.0-albedo[i]);


	    }

	}

	albedoCount = albedoCount +1;
    }

void World::calcSurfaceAlbedo(Body* star, int iLatitude, double meanZenith)
    {
    /*
     * Written 01/05/2017 by dh4gan
     * Calculates the surface albedo
     *
     */

    double aIce,aOcean;

    // Calculate the CO2 saturation vapour pressure

    double psat = 0.0;
    if(T[iLatitude]< TsatSolid)
	{
	psat = 6.760956 - 1284.07/(T[iLatitude] - 4.718) + 1.256E-4*(T[iLatitude] - 143.15);
	}

    else
	{
	psat = 3.128082 - 867.2124/T[iLatitude] + 1.865612E-2*T[iLatitude] - 7.248820E-5*T[iLatitude]*T[iLatitude] + 9.3E-8*T[iLatitude]*T[iLatitude]*T[iLatitude];
	}
    psat = pow(10.0, psat);
    psat = 1.013*psat;


    // If the local CO2 pressure exceeds this saturation pressure, then CO2 freezes out

    if(CO2pressure[iLatitude]>psat)
	{

	// Albedo = CO2 ice albedo
	surfaceAlbedo[iLatitude] = aCO2Ice;

	}

    else
	{
	// Calculate ice albedo
	aIce = aIceVisible*star->getfVisible() + aIceIR*star->getfIR();

	// Calculate ocean albedo (WK97)
	double aCloud = -0.078 + 0.65*meanZenith;
	aOcean = 0.8; // TODO - get true ocean albedo (email Jacob)

	surfaceAlbedo[iLatitude] = (aLand*landFraction + aOcean*oceanFraction)*(1.0-iceFraction[iLatitude]) +
	    aIce*iceFraction[iLatitude];

	//printf("AS: %e %e %e %e \n", surfaceAlbedo[iLatitude], aLand, aOcean, aIce);
	}

    }
void World::calcAlbedo(Body* star, int iLatitude, double meanZenith)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Simple calculation of the albedo as a function of latitude according to Temperature
     *
     */

    if (CSCycleOn)
	{

	// If CS Cycle on, use fitting function from Haqq-Misra et al (2016)
	calcSurfaceAlbedo(star, iLatitude,meanZenith);

	vector<double> coeff = star->getAlbedoCoefficients(T[iLatitude]);

	double as = surfaceAlbedo[iLatitude];
	double phi = log10(CO2pressure[iLatitude]);
	double logT = log10(T[iLatitude]);
	double mu = meanZenith;

	albedo[iLatitude] =  coeff[0]*mu*mu*mu + coeff[1]*mu*mu*as +
		coeff[2]*mu*mu*logT + coeff[3]*mu*mu*phi + coeff[4]*mu*mu +
		coeff[5]*mu*as*as + coeff[6]*mu*as*logT + coeff[7]*mu*as*phi +
		coeff[8]*mu*as + coeff[9]*mu*logT*logT + coeff[10]*mu*logT*phi +
		coeff[11]*mu*logT + coeff[12]*mu*phi*phi + coeff[13]*mu*phi +
		coeff[14]*mu + coeff[15]*as*as*as + coeff[16]*as*as*logT +
		coeff[17]*as*as*phi + coeff[18]*as*as +coeff[19]*as*logT*logT +
		coeff[20]*as*logT*phi + coeff[21]*as*logT + coeff[22]*as*phi*phi +
		coeff[23]*as*phi + coeff[24]*as + coeff[25]*logT*logT*logT +
		coeff[26]*logT*logT*phi + coeff[27]*logT*logT + coeff[28]*logT*phi*phi +
		coeff[29]*logT*phi + coeff[30]*logT + coeff[31]*phi*phi*phi +
		coeff[32]*phi*phi + coeff[33]*phi + coeff[34];

	meanAlbedo[iLatitude] = meanAlbedo[iLatitude] + albedo[iLatitude];
	}

    else
	{
	albedo[iLatitude] = 0.525 - 0.245 * tanh((T[iLatitude] - freeze + 5.0) / 5.0);
	meanAlbedo[iLatitude] = albedo[iLatitude];
	albedoCount = 1;
	}




    }

void World::calcHeatCapacity(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Calculation of Heat Capacity of the atmosphere as a function of the local temperature
     *
     */
    double C_ice = 0.0;
    double C_land = 5.25e9; //differs by factor 1000 to W&K...
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


void World::calcCO2Rates(int iLatitude)
    {
    /*
     * Written 01/05/2017 by dh4gan
     * Calculates the rate of change of CO2 partial pressure using regime defined
     * in Haqq-Misra et al (2016)
     *
     */

    // Calculate Weathering Rate

    landWeathering[iLatitude] = pow(CO2pressure[iLatitude]/CO2Earth, betaCO2)*(exp(kactive*(T[iLatitude]-288.0)));
    landWeathering[iLatitude] = outgassingRate*landWeathering[iLatitude]*(1.0 + krun*(T[iLatitude]-288.0))/year;

    oceanWeathering[iLatitude] = W0*pow(CO2pressure[iLatitude]/CO2Earth, gammaCO2);


    CO2dot[iLatitude] = outgassingRate - landWeathering[iLatitude] - oceanWeathering[iLatitude];
    //printf("CO2 dot: %e %e %e %e %e %e %e %e \n", CO2dot[iLatitude], T[iLatitude]-288.0,betaCO2,kactive,krun, outgassingRate, landWeathering[iLatitude],oceanWeathering[iLatitude]);

	   
    }





void World::calcCooling(int iLatitude)
    {
    /*
     * Written 9/1/14 by dh4gan
     * Edited by dh4gan 01/05/2017
     * Calculates the Infrared Cooling as a function of optical depth and temperature
     *
     */

	if(CSCycleOn)
		{
		double phi = log10(CO2pressure[iLatitude]);
		double logT = log10(T[iLatitude]); //

		// Cooling function in (W m^-2)
		infrared[iLatitude] = ircoeff[0]*logT*logT*logT*logT +
				    ircoeff[1]*logT*logT*logT*phi +
				    ircoeff[2]*logT*logT*logT +
				    ircoeff[3]*logT*logT*phi*phi +
				    ircoeff[4]*logT*logT*phi +
				    ircoeff[5]*logT*logT +
				    ircoeff[6]*logT*phi*phi*phi +
				    ircoeff[7]*logT*phi*phi +
				    ircoeff[8]*logT*phi +
				    ircoeff[9]*logT +
				    ircoeff[10]*phi*phi*phi*phi +
				    ircoeff[11]*phi*phi*phi +
				    ircoeff[12]*phi*phi +
				    ircoeff[13]*phi +
				    ircoeff[14];

		// convert to erg s^-1 cm^-2
		infrared[iLatitude] = infrared[iLatitude]*1000.0;


		//printf("%i %f %f %f %f \n", iLatitude, CO2pressure[iLatitude],CO2Earth,logT, infrared[iLatitude]);
		}
	else{
		infrared[iLatitude] = sigma_SB*pow(T[iLatitude],4)/(1.0+0.75*tau[iLatitude]); //erg s^-1 cm^-2
	}
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

	Q[iLatitude] = absorbedInsol[iLatitude] + tidal[iLatitude]- infrared[iLatitude];
	//printf("Q: %e %e %e \n", absorbedInsol[iLatitude], tidal[iLatitude], infrared[iLatitude]);
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

	    //select diffusion coeff for following calculations
    	    double diffusionvalue;
    	    if (CSCycleOn){diffusionvalue = diffusion[i];}
    	    else {diffusionvalue = diffusion0;}

	    if (i == nPoints)
		{
		Dplus = diffusionvalue * (1.0 - x[i] * x[i]);
		}
	    else
		{
		Dplus = diffusionvalue * 0.5
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
		cout << timestep[i] << "   " << T[i] << "   " << C[i] << diffusion[i] << endl;
	}

	cout << dtLEBM <<  endl;
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
	
	    //select diffusion coeff for following calculations
    	    double diffusionvalue;
    	    if (CSCycleOn){diffusionvalue = diffusion[i];}
    	    else {diffusionvalue = diffusion0;}

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
		Dminus = diffusionvalue * (1.0 - x[i] * x[i]);
		}
	    else
		{
		Tminus1 = T_old[i - 1];
		Dminus = 0.5 * diffusionvalue
			* ((1.0 - x[i - 1] * x[i - 1]) + (1.0 - x[i] * x[i]));
		}

	    // Now the maximum value

	    if (i == nPoints)
		{
		Tplus1 = T_old[i];
		Dplus = diffusionvalue * (1.0 - x[i] * x[i]);
		}
	    else
		{
		Tplus1 = T_old[i + 1];
		Dplus = 0.5 * diffusionvalue
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


	    // Update the CO2 pressure (rates /second)

	    if(CSCycleOn)
		{
		CO2pressure[i] = CO2pressure[i] + CO2dot[i]*dtLEBM;
		printf("CO2 update: %e %e %e \n", CO2pressure[i], CO2dot[i], dtLEBM);
		}


	    // Set up diffusion constant
	    diffusion[i] = 5.8e2*(CO2pressure[i]/CO2Earth)*pow(rotationPeriod, 2.);
		//erg cm^-2 s^-1 K^-1
		//W&K expression


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
	double &meanhab, double &meanTidal, double &meanCO2p, double &meanDiffusion)
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
    meanCO2p = 0.0;
    meanDiffusion = 0.0;

    for (int i=0; i< nPoints1; i++)
	{
	meanT = meanT + T[i]*0.5*deltax[i];
	meanQ = meanQ + Q[i]*0.5*deltax[i];
	meanA = meanA + meanAlbedo[i]*0.5*deltax[i];
	meanIR = meanIR + infrared[i]*0.5*deltax[i];
	meanS = meanS + insol[i]*0.5*deltax[i];
	meanhab = meanhab + hab[i]*0.5*deltax[i];
	meanTidal = meanTidal + tidal[i]*0.5*deltax[i];
	meanCO2p = meanCO2p + CO2pressure[i]*0.5*deltax[i];
	meanDiffusion = meanDiffusion + diffusion[i]*0.5*deltax[i];

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

    double minT, maxT, meanT, meanQ, meanA, meanIR, meanS, meanhab, meanTidal,meanCO2p, meanDiffusion;
    string formatString;

    // Get meanAlbedo by dividing by number of objects contributing

    for (int i=0; i<nPoints1; i++)
	{
	meanAlbedo[i] = meanAlbedo[i]/float(albedoCount);
	}

    // Firstly, write line to log file
    calcLEBMMeans(minT, maxT, meanT, meanQ, meanA, meanIR,meanS, meanhab, meanTidal, meanCO2p, meanDiffusion);

    // Also include orbital data here

    formatString = "+%.8E  ";

    for (int icol=0;icol < 19; icol++)
	{
	formatString = formatString + "%.8E  ";
	}
    formatString = formatString + "\n";

    fprintf(logFile, formatString.c_str(),
	    tSnap, minT, maxT, meanT,meanQ,meanA,meanIR,meanS,meanhab, meanTidal,
	    semiMajorAxis, eccentricity, inclination,
	    argumentPeriapsis, longitudeAscendingNode, meanAnomaly,
	    obliquity*radToDeg,precession*radToDeg,meanCO2p, meanDiffusion);

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
	fprintf(snapshotFile, "%+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E %+.6E \n",
		x[i], lat[i], T[i], C[i], Q[i], infrared[i],
		meanAlbedo[i], insol[i], tau[i],iceFraction[i], hab[i], tidal[i], CO2pressure[i], diffusion[i]);
	}
    fclose(snapshotFile);
    }

    }

int World::getRestartParameters()
    {

    /*
     * written by dh4gan 13/5/14
     * Reads temperature data and obliquity/precession from last file entry for restart purposes
     *
     */

    // Read log file to determine how many snapshots there have been

    string logFileName = getName()+".log";
    string latFileName = getName()+".lat";
    string line;
    int nSnap = 0;
    double dtmax = 1.0e30;


    int i;
    double d;

    // Obtain temperature data from last line of .lat file
    ifstream latfile(latFileName.c_str());

    if(latfile)
	{
	while(getline(latfile,line))
	    {
	    i=0;
	    istringstream iss(line);
	     while (iss >> d) {
	        T[i] =d;
	        i++;
	     }
	    }
	}
    else
       	{
       	cout << "Error: File "<< latFileName << "not found" << endl;
       	return -1;
       	}


    // Read last line of .log file to obtain obliquity / precession data

    ifstream lfile(logFileName.c_str());

    if(lfile)
   	{

   	while(getline(lfile,line))
   	    {
   	    nSnap++;
   	    }

   	cout << "Number of Snapshots is " << nSnap << endl;

   	// Extract data from last line

   	istringstream iss(line);

   	// Skip the other data
   	for (int icol=0; icol<<15; icol++)
   	    {
   	    iss>>d;
   	    }
   	iss >> obliquity;
   	iss >> precession;

   	printf("Obliquity read as %f degrees \n",obliquity*radToDeg);
   	printf("Precession read as %f degrees \n",precession*radToDeg);


   	}
       else
   	{
   	cout << "Error: File "<< logFileName << "not found" << endl;
   	return -1;
   	}

    lfile.close();

    // Recalculate EBM data model

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
	    //calcAlbedo(i);
	    calcCooling(i);
	    calcNetHeating(i);
	    calcHabitability(i,freeze,boil);
	    }
	}
   calcLEBMTimestep(dtmax);

   return nSnap;
    }



