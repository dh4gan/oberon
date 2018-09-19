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
#include <string.h>
#include <map>
#include "System.h"
#include "Star.h"

#include <fstream>
#include <sstream>
//using namespace std;

// TODO - rewrite parFile class using STL map objects

const string stringType = "string";
const string intType = "int";
const string doubleType = "double";
const string vectorIntType = "vectorInt";
const string vectorDoubleType = "vectorDouble";
const string vectorStringType = "vectorString";

class parFile {
public:
	parFile( );
	parFile(string name);
    

    
    
	string NBodyFile;
	string parFileName;
	string SystemName;
	string fileType;

    
    
    
    // Map objects to store variable data as it is read
    
    std::map < string, string > stringVariables;
    std::map < string, double > doubleVariables;
    std::map < string, int > intVariables;
    
    
    std::map < string, vector<int> > vectorIntVariables;
    std::map < string, vector<double> > vectorDoubleVariables;
    std::map < string, vector<string> > vectorStringVariables;
    
    // This map stores which map each variable lives in
    std::map< string, string > variableLocations;
    
    
	bool restart;
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

	Vector3D getBodyPosition(int index);
	Vector3D getBodyVelocity(int index);

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
    
    // TODO - write initialiseVectors and initialiseBooleans
    void initialiseVectors(int nBodies);
    void initialiseBooleans();
    
    string getStringVariable(string &par){return stringVariables[par];}
    string getStringVariable(string &par, int bodyIndex){return vectorStringVariables[par][bodyIndex];}
    
    int getIntVariable(string &par){return intVariables[par];}
    int getIntVariable(string &par, int &bodyIndex){return vectorIntVariables[par][bodyIndex];}
    double getDoubleVariable(string &par){return doubleVariables[par];}
    double getDoubleVariable(string &par, int &bodyIndex){return vectorDoubleVariables[par][bodyIndex];}
    
	void setupRestartPositions();


};


#endif /* PARFILE_H_ */
