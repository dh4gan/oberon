/*
 * parFile.h
 *
 *  Created on: Sep 23, 2013
 *  Majorly revised: Sep 19, 2018
 *
 *      Author: dh4gan
 *
 * This object reads in the input parameters for the simulation
 * Parameters stored in <map> objects
 *
 */

#ifndef PARFILE_H_
#define PARFILE_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include "Vector3D.h"

#include <fstream>
#include <sstream>
using namespace std;

const string stringType = "string";
const string intType = "int";
const string doubleType = "double";
const string vectorIntType = "vectorInt";
const string vectorDoubleType = "vectorDouble";
const string vectorStringType = "vectorString";


// All single string variables to read in
const string stringV[] = {"ParType", "NBodyOutput", "SystemName","Restart","ObliquityEvolution", "CarbonateSilicateCycle", "TidalHeating","PlanetaryIllumination","FullOutput"};

// Of these, which are boolean?
 const string boolV[] = {"Restart","ObliquityEvolution", "IceMeltingOn","CarbonateSilicateCycle", "TidalHeating","PlanetaryIllumination","FullOutput"};

// Scalar ints to be read in
const string intV[] = {"NGridPoints","Number_Bodies"};

// Scalar doubles
const string doubleV[] = {"TotalMass","SnapshotTime","MaximumTime"};

// Vector (string) variables
const string vectorStringV[] = {"BodyName", "BodyType", "IceMeltingOn","SpectralType"};

// Vector (int) variables
const string vectorIntV[] = {"OrbitCentre"};

// Vector (double) variables
const string vectorDoubleV[] = {"Mass", "Radius", "Position", "XPosition", "YPosition", "ZPosition", "Velocity", "XVelocity", "YVelocity", "ZVelocity", "SemiMajorAxis", "Eccentricity", "Inclination", "LongAscend", "Periapsis", "MeanAnomaly", "RotationPeriod", "Obliquity", "WinterSolstice", "OceanFraction", "InitialTemperature", "Albedo","Luminosity", "OutgassingRate", "BetaCO2", "GammaCO2","SeafloorWeatheringRate",};


const string degreeV[] = {"Obliquity","WinterSolstice"}; // Variables with units of degrees - these are eventually converted to radians

// Decant these into STL vectors for easier passing between functions
const vector<const string> stringVar(stringV, stringV+sizeof(stringV)/sizeof(*stringV));
const vector<const string> boolVar(boolV, boolV+sizeof(boolV)/sizeof(*boolV));
const vector<const string> intVar(intV, intV+sizeof(intV)/sizeof(*intV));
const vector<const string> doubleVar(doubleV, doubleV+sizeof(doubleV)/sizeof(*doubleV));

const vector<const string> vectorStringVar(vectorStringV, vectorStringV+sizeof(vectorStringV)/sizeof(*vectorStringV));
const vector<const string> vectorIntVar(vectorIntV, vectorIntV+sizeof(vectorIntV)/sizeof(*vectorIntV));
const vector<const string> vectorDoubleVar(vectorDoubleV, vectorDoubleV+sizeof(vectorDoubleV)/sizeof(*vectorDoubleV));

const vector<const string> degreeVar(degreeV, degreeV+sizeof(degreeV)/sizeof(*degreeV));

class parFile {
public:
	parFile( );
	parFile(string name);
    
    virtual ~parFile();
    
    void readFile();
    void readFile(string &filename);
    
    void readVariable(string &par, istringstream &iss, int &bodyIndex);
    void readIntVariable(string &par, istringstream &iss);
    void readStringVariable(string &par, istringstream &iss);
    void readDoubleVariable(string &par, istringstream &iss);
    
    void readVectorIntVariable(string &par, istringstream &iss, int &bodyIndex);
    void readVectorDoubleVariable(string &par, istringstream &iss, int &bodyIndex);
    void readVectorStringVariable(string &par, istringstream &iss, int &bodyIndex);
    void read3DVector(string &par, istringstream &iss, int &bodyIndex);
    
    string getStringVariable(const string &par){return stringVariables[par];}
    string getStringVariable(const string &par, int &bodyIndex){return vectorStringVariables[par][bodyIndex];}
    
    int getIntVariable(const string &par){return intVariables[par];}
    int getIntVariable(const string &par, int &bodyIndex){return vectorIntVariables[par][bodyIndex];}
    double getDoubleVariable(const string &par){return doubleVariables[par];}
    double getDoubleVariable(const string &par, int &bodyIndex){return vectorDoubleVariables[par][bodyIndex];}
    
    bool getBoolVariable(const string &par){return boolVariables[par];}
    
    void setVariableType(const vector<const string> &variables, const string &type);
    void setVariableLocations();
    
    Vector3D getBodyPosition(int index);
    Vector3D getBodyVelocity(int index);
    
    void initialiseVectors(int nBodies);
    void convertToRadians(int nBodies);
    
    void initialiseBoolean(const string &par);
    void initialiseAllBooleans();
    
    void setupRestartPositions();
    
    void checkParameters(); // TODO - write checkSetup routine
    void displayParameters();
    void reportError(const string &par, double &value);
    void reportError(const string &par, int &value);
    void reportError(const string &par, string &value);
    
private:
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
    
	
    
};


#endif /* PARFILE_H_ */
