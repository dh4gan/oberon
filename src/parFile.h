/*
 * parFile.h
 *
 *  Created on: Sep 23, 2013
 *      Author: davidharvey
 */

#ifndef PARFILE_H_
#define PARFILE_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
//#include "System.h"
//#include "Star.h"
#include "Vector3D.h"

#include <fstream>
#include <sstream>
using namespace std;

// TODO - rewrite parFile class using STL map objects

const string stringType = "string";
const string intType = "int";
const string doubleType = "double";
const string vectorIntType = "vectorInt";
const string vectorDoubleType = "vectorDouble";
const string vectorStringType = "vectorString";


// All single string variables to read in
const string stringVar[] = {"ParType", "NBodyOutput", "SystemName","Restart","ObliquityEvolution", "CarbonateSilicateCycle", "TidalHeating","PlanetaryIllumination","FullOutput"};

// Of these, which are boolean?
 const string boolVar[] = {"Restart","ObliquityEvolution", "IceMeltingOn","CarbonateSilicateCycle", "TidalHeating","PlanetaryIllumination","FullOutput"};

// Scalar ints
const string intVar[] = {"NGridPoints","Number_Bodies"};

// Scalar doubles
const string doubleVar[] = {"TotalMass","SnapshotTime","MaximumTime"};

// Vector (string) variables
const string vectorStringVar[] = {"BodyName", "BodyType", "IceMeltingOn"};

// Vector (int) variables
const string vectorIntVar[] = {"orbitCentre"};

// Vector (double) variables
const string vectorDoubleVar[] = {"Mass", "Radius", "Position", "XPosition", "YPosition", "ZPosition", "Velocity", "XVelocity", "YVelocity", "ZVelocity", "SemiMajorAxis", "Eccentricity", "Inclination", "LongAscend", "Periapsis", "MeanAnomaly", "RotationPeriod", "Obliquity", "WinterSolstice", "OceanFraction", "InitialTemperature", "Luminosity"};

class parFile {
public:
	parFile( );
	parFile(string name);
    
	string parFileName;
	
    // Map objects to store variable data as it is read
    
    std::map < string, string > stringVariables;
    std::map < string, double > doubleVariables;
    std::map < string, int > intVariables;
    std::map < string, bool > boolVariables;
    
    
    std::map < string, vector<int> > vectorIntVariables;
    std::map < string, vector<double> > vectorDoubleVariables;
    std::map < string, vector<string> > vectorStringVariables;
    
    // This map stores which map each variable lives in
    std::map< string, string > variableLocations;
    
	/*bool restart;
	bool illuminationOn;
	bool obliquityOn;
	bool tidalHeatingOn;
    bool CSCycleOn;
	bool fullOutput;

	int number_bodies;
	int nPoints;
	int snapshotNumber;
	double snapshotTime;
	double maximumTime;
	double systemTime;
	double totalMass;

	vector<string> parameters;
	vector<string> BodyNames;
	vector<string> BodyTypes;

	vector<double> Mass;
	vector<double> Radius;

	vector<double> x_position;
	vector<double> y_position;
	vector<double> z_position;

	vector<double> x_velocity;
	vector<double> y_velocity;
	vector<double> z_velocity;

	vector<double> semiMajorAxis;
	vector<double> eccentricity;
	vector<double> inclination;
	vector<double> longAscend;
	vector<double> Periapsis;
	vector<double> meanAnomaly;

	vector<string> spectralType;
	vector<double> luminosity;
	vector<double> albedo;

	vector<double> rotationPeriod;
	vector<double> obliquity;
	vector<double> precession;
	vector<double> oceanFraction;
	vector<double> initialTemperature;
	vector<double> outgassingRate;
	vector<double> betaCO2;
	vector<double> gammaCO2;
	vector<double> seafloorWeathering;
	vector<bool> activateMelt;

	vector<int> orbitCentre;
*/
	Vector3D getBodyPosition(int index);
	Vector3D getBodyVelocity(int index);
    
    void readFile();
    void readFile(string &filename);
    
	int readParFile();
	int readParFile(string fileName);

	void readPosFile();
	void readOrbFile();
	int parType();
	int parType(string fileName);

    
    void setVariableLocations();
      
    void readVariable(string &par, istringstream &iss, int &bodyIndex);
    void readIntVariable(string &par, istringstream &iss);
    void readStringVariable(string &par, istringstream &iss);
    void readDoubleVariable(string &par, istringstream &iss);
    
    void read3DVector(string &par, istringstream &iss, int &bodyIndex);
    
    void readVectorVariable(string &par, istringstream &iss);
    
    void readVectorIntVariable(string &par, istringstream &iss, int &bodyIndex);
    void readVectorDoubleVariable(string &par, istringstream &iss, int &bodyIndex);
    void readVectorStringVariable(string &par, istringstream &iss, int &bodyIndex);
    
    void initialiseVectors(int nBodies);
    void convertToRadians(int nBodies);
    
    bool initialiseBoolean(string &par);
    
    void initialiseAllBooleans();
    
    
    string getStringVariable(const string &par){return stringVariables[par];}
    string getStringVariable(const string &par, int &bodyIndex){return vectorStringVariables[par][bodyIndex];}
    
    int getIntVariable(const string &par){return intVariables[par];}
    int getIntVariable(const string &par, int &bodyIndex){return vectorIntVariables[par][bodyIndex];}
    double getDoubleVariable(const string &par){return doubleVariables[par];}
    double getDoubleVariable(const string &par, int &bodyIndex){return vectorDoubleVariables[par][bodyIndex];}
    
    bool getBoolVariable(const string &par){return boolVariables[par];}
   // bool getBoolVariable(string &par, index &bodyIndex){return boolVariables[par][bodyIndex];}

    
	void setupRestartPositions();
    
    void checkSetup(); // TODO - write checkSetup routine
    
    void printParameterOptions();// TODO - write Parameter options to screen (file?)


    // TODO - create Body Objects from parFile object
    
};


#endif /* PARFILE_H_ */
