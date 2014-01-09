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
#include "System.h"
#include "Star.h"

#include <fstream>
#include <sstream>
using namespace std;

class parFile {
public:
	parFile( );
	parFile(string name);

	string outputPrefix;
	string parFileName;
	string SystemName;
	string fileType;
	int number_bodies;
	double outputFrequency;
	double maximumTime;
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
	vector<double> luminosity;
	vector<double> rotationPeriod;
	vector<double> obliquity;
	vector<double> winterSolstice;
	vector<double> oceanFraction;
	vector<double> initialTemperature;



	Vector3D getBodyPosition(int index);
	Vector3D getBodyVelocity(int index);

	int readParFile();
	int readParFile(string fileName);

	void readPosFile();
	void readOrbFile();
	int parType();
	int parType(string fileName);


};


#endif /* PARFILE_H_ */
