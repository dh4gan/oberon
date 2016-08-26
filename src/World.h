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

    World(string namestring, double m, double rad,Vector3D pos, Vector3D vel, int n, double obliq, double rot, double prec,
	    double ocean, double T0, bool melt, bool start, bool tide, bool obevol, bool CScycle);

    World(string namestring, double m, double rad,
	    double semimaj, double ecc, double inc, double longascend,
	    double argper, double meananom, double G, double totalMass,
	    int n, double obliq, double rot, double prec,
	    double ocean, double T0, bool melt, bool start, bool tide, bool obevol, bool CScycle);

    virtual ~World();

    /* Accessors */

    void setInsolationZero() {insol.assign(nPoints1,0.0);}
    void setTemperature(vector<double> temp);

    //void setHostBody(Body* bod){hostBody = bod;}

    double getLEBMTimestep(){return dtLEBM;}

    double getObliquity(){return obliquity;}

    // Calculation Methods

    void initialiseLEBM();
    void updateLEBM(vector<Body*> bodies, double &G, double &totmass, vector<double>eclipsefrac,bool &planetaryIllumination);
    void updateLEBM(vector<Body*> bodies, double &G, double &totmass, vector<double> eclipsefrac, double &dtmax, bool &planetaryIllumination);

    void calcObliquity(vector<Body*>bodies, double G, double totmass);
    void calcInsolation(Body* star, double &eclipsefrac);
    void calcAlbedo(int iLatitude);
    void calcHeatCapacity(int iLatitude);
    void calcIce(int iLatitude);
    void calcOpticalDepth(int iLatitude);
    void calcCooling(int iLatitude);
    void calcTidalHeating(int iLatitude);
    void calcCO2pressure(int iLatitude);
    void calcNetHeating(int iLatitude);
    void calcHabitability(int iLatitude,double &minT, double &maxT);
    void calcLEBMTimestep(double &dtmax);

    void integrate();

    void calcLuminosity();

    double getLuminosity(){return luminosity;}

    // Output Methods

    void outputLEBMData(int &snapshotNumber, double &tSnap, bool fullOutput);
    void initialiseOutputVariables(bool restart);
    void calcLEBMMeans(double &minT, double &maxT, double &meanT, double &meanQ, double &meanA, double &meanIR, double &meanS,
    	double &meanhab, double &meanTidal, double &meanCO2p, double &meanDiffusion);

    int findRestartTemperature();

    // Standard cloning method
    virtual World* Clone()
	{
	return new World(*this);
	}

protected:
    double rotationPeriod;
    double obliquity;
    double ellipticity;
    double precession;
    double oceanFraction;
    double landFraction;
    double initialTemperature;
    double diffusion0;

    double luminosity;

    double rho_moon;
    double rigid;
    double Qtidal;


    bool CSCycleOn; //Giblin 10/7/15.


    int nPoints,nPoints1;

    double nFloat;
    double dtLEBM;

    vector<double> lat, x,coslat,tanlat, deltax;
    vector<double> T, T_old, tau;
    vector<double> iceFraction, C, hab;
    vector<double> infrared, Q, albedo,tidal, insol, CO2pressure, diffusion;

    FILE *logFile, *snapshotFile, *latFile;
    bool activateMelt;
    bool obliquityEvolutionOn;
    bool restart;

    bool tidalHeatingOn;
    vector<bool> melting;
    vector<double> meltTime;

    };

#endif /* WORLD_H */
