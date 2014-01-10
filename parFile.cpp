/*
 * parFile.cpp
 *
 *  Created on: Sep 23, 2013
 *
 *      Author: davidharvey
 */

#include "parFile.h"
#include <iostream> // Included for debug lines only
#include <math.h>
#include <string>

parFile::parFile()
    {

    }

parFile::parFile(string name)
    {
    parFileName = name;

    }

Vector3D parFile::getBodyPosition(int index)
    {
    /*
     * Author: dh4gan
     * Extracts a Vector3D object containing the position of body index
     *
     */

    Vector3D position;

    position.elements[0] = x_position[index];
    position.elements[1] = y_position[index];
    position.elements[2] = z_position[index];

    return position;

    }

Vector3D parFile::getBodyVelocity(int index)
    {
    /*
     * Author: dh4gan
     * Extracts a Vector3D object containing the position of body index
     *
     */

    Vector3D velocity;

    velocity.elements[0] = x_velocity[index];
    velocity.elements[1] = y_velocity[index];
    velocity.elements[2] = z_velocity[index];

    return velocity;

    }

void parFile::readPosFile()
    {
    /*
     * Read par files from input and set up the body and parse back
     * body
     *
     */
    // First read in the user data file
    // Format must be,name of body followed by its parameters

    string par;
    string line;
    string BodyType;

    char inputfile[100];

    double val_i;
    double val_k;
    double val_j;

    strcpy(inputfile, parFileName.c_str());

    ifstream myfile(inputfile);

    // Then loop through each line using getline and then
    //assign to vectors

    // Handy defaults for incomplete files

    luminosity.push_back(0.0);
    rotationPeriod.push_back(0.0);
    obliquity.push_back(0.0);
    winterSolstice.push_back(0.0);
    oceanFraction.push_back(0.0);
    initialTemperature.push_back(0.0);

    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;

	cout << par << endl;

	if(par=="OutputPrefix")
		    {
		    iss >> outputPrefix;
		    }
		if(par=="OutputFrequency")
		    {
		    iss >> outputFrequency;
		    }

		if(par =="MaximumTime")
		    {
		    iss >> maximumTime;
		    }

	if (par == "SystemName")
	    {
	    iss >> SystemName;
	    }
	if (par == "Number_Bodies")
	    {
	    iss >> number_bodies;

	    }
	else if(par =="NGridPoints")
	    {
	    iss >> nPoints;
	    }
	else if (par == "BodyType")
	    {
	    iss >> BodyType;
	    BodyTypes.push_back(BodyType);
	    }
	else if (par == "Mass")
	    {
	    iss >> val_i;
	    Mass.push_back(val_i);
	    totalMass += val_i;
	    }
	else if (par == "Radius")
	    {
	    iss >> val_i;
	    Radius.push_back(val_i);
	    }
	else if (par == "Position")
	    {
	    iss >> val_i >> val_j >> val_k;
	    x_position.push_back(val_i);
	    y_position.push_back(val_j);
	    z_position.push_back(val_k);
	    }
	else if (par == "Velocity")
	    {
	    iss >> val_i >> val_j >> val_k;
	    x_velocity.push_back(val_i);
	    y_velocity.push_back(val_j);
	    z_velocity.push_back(val_k);
	    }

	else if (par == "Luminosity")
	    {
	    iss >> val_i;
	    luminosity.pop_back(); // Discard default
	    luminosity.push_back(val_i);
	    }

	else if (par == "RotationPeriod")
	    {
	    iss >> val_i;
	    rotationPeriod.pop_back(); // Discard default
	    rotationPeriod.push_back(val_i);
	    }

	else if (par == "Obliquity")
	    {
	    iss >> val_i;
	    obliquity.pop_back(); // Discard default
	    obliquity.push_back(val_i);
	    }

	else if (par == "WinterSolstice")
	    {
	    iss >> val_i;
	    winterSolstice.pop_back(); // Discard default
	    winterSolstice.push_back(val_i);
	    }

	else if (par == "OceanFraction")
	    {
	    iss >> val_i;
	    oceanFraction.pop_back(); // Discard default
	    oceanFraction.push_back(val_i);
	    }
	else if (par == "InitialTemperature")
	    {
	    iss >> val_i;
	    initialTemperature.pop_back(); // Discard default
	    initialTemperature.push_back(val_i);
	    }

	}
    myfile.close();
    }
void parFile::readOrbFile()
    {
    /*
     * Read par files from input and set up the body and parse back
     * body
     *
     */
    // First read in the user data file
    // Format must be,name of body followed by its parameters

    string par;
    string line;
    string BodyType, BodyName;

    char inputfile[100];

    double val_i;

    strcpy(inputfile, parFileName.c_str());

    ifstream myfile(inputfile);

    // Set some handy defaults here if the input file is incomplete

    longAscend.push_back(0.0);
    inclination.push_back(0.0);
    Periapsis.push_back(0.0);
    meanAnomaly.push_back(0.0);

    luminosity.push_back(0.0);
    rotationPeriod.push_back(0.0);
    obliquity.push_back(0.0);
    winterSolstice.push_back(0.0);
    oceanFraction.push_back(0.0);
    initialTemperature.push_back(0.0);

    // Then loop through each line using getline and then
    //assign to vectors
    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;

	if(par=="OutputPrefix")
	    {
	    iss >> outputPrefix;
	    }
	if(par=="OutputFrequency")
	    {
	    iss >> outputFrequency;
	    }

	if(par =="MaximumTime")
	    {
	    iss >> maximumTime;
	    }

	if (par == "SystemName")
	    {
	    iss >> SystemName;
	    }
	if (par == "Number_Bodies")
	    {
	    iss >> number_bodies;

	    }

	else if(par =="NGridPoints")
		    {
		    iss >> nPoints;
		    }

	else if (par == "BodyName")
	    {
	    iss >> BodyName;
	    BodyNames.push_back(BodyName);
	    }

	else if (par == "BodyType")
	    {
	    iss >> BodyType;
	    BodyTypes.push_back(BodyType);
	    }
	else if (par == "Mass")
	    {
	    iss >> val_i;
	    Mass.push_back(val_i);
	    totalMass += val_i;
	    }
	else if (par == "Radius")
	    {
	    iss >> val_i;
	    Radius.push_back(val_i);
	    }
	else if (par == "SemiMajorAxis")
	    {
	    iss >> val_i;
	    semiMajorAxis.push_back(val_i);
	    }
	else if (par == "Eccentricity")
	    {
	    iss >> val_i;
	    eccentricity.push_back(val_i);
	    }
	else if (par == "Inclination")
	    {
	    iss >> val_i;
	    inclination.pop_back(); // Discard default
	    inclination.push_back(val_i);
	    }
	else if (par == "LongAscend")
	    {
	    iss >> val_i;
	    longAscend.pop_back(); // Discard default
	    longAscend.push_back(val_i);
	    }

	else if (par == "Periapsis")
	    {
	    iss >> val_i;
	    Periapsis.pop_back(); // Discard default
	    Periapsis.push_back(val_i);
	    }
	else if (par == "MeanAnomaly")
	    {
	    iss >> val_i;
	    meanAnomaly.pop_back(); // Discard default
	    meanAnomaly.push_back(val_i);
	    }

	else if (par == "Luminosity")
	    {
	    iss >> val_i;
	    luminosity.pop_back(); // Discard default
	    luminosity.push_back(val_i);
	    }

	else if (par == "RotationPeriod")
	    {
	    iss >> val_i;
	    rotationPeriod.pop_back(); // Discard default
	    rotationPeriod.push_back(val_i);
	    }

	else if (par == "Obliquity")
	    {
	    iss >> val_i;
	    obliquity.pop_back(); // Discard default
	    obliquity.push_back(val_i);
	    }

	else if (par == "WinterSolstice")
	    {
	    iss >> val_i;
	    winterSolstice.pop_back(); // Discard default
	    winterSolstice.push_back(val_i);
	    }

	else if (par == "OceanFraction")
	    {
	    iss >> val_i;
	    oceanFraction.pop_back(); // Discard default
	    oceanFraction.push_back(val_i);
	    }
	else if (par == "InitialTemperature")
	    {
	    iss >> val_i;
	    initialTemperature.pop_back(); // Discard default
	    initialTemperature.push_back(val_i);
	    }

	}
    myfile.close();
    }
int parFile::parType()
    {
    /*
     * PURPOSE : TO DETERMINE THE TYPE OF PARAMETER FILE
     *           BEING USED. REQUIRES THE INPUT FILE TO HAVE
     *           A HEADER STATING WHICH FILE TO BE USED
     *
     * RETURN : INTEGER :
     * 				0 = POSITION INFORMATION
     * 				1 = ORBITAL PARAMETERS
     *
     *
     */
    char inputfile[100];
    string parType;
    string line;
    string par;

    printf("Reading in the Users body data \n ");
    cout << "What is the input file? " << endl;

    getline(cin, parFileName);

    cout << "Reading input from file " << parFileName << endl;

    strcpy(inputfile, parFileName.c_str());
    ifstream myfile(inputfile);

    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;
	if (par == "ParType")
	    {
	    iss >> parType;
	    }
	}
    if (parType == "Positional")
	{
	fileType = parType;
	return 0;
	}
    else if (parType == "Orbital")
	{
	fileType = parType;
	return 1;
	}
    else
	{
	return 2;
	}
    }

int parFile::parType(string fileName)
    {
    /*
     * PURPOSE : TO DETERMINE THE TYPE OF PARAMETER FILE
     *           BEING USED. ACCEPTS INPUT FILE AS ARGUMENT.
     *           REQUIRES THE INPUT FILE TO HAVE
     *           A HEADER STATING WHICH FILE TO BE USED
     *
     * RETURN : INTEGER :
     * 				0 = POSITION INFORMATION
     * 				1 = ORBITAL PARAMETERS
     *
     *
     */
    char inputfile[100];
    string parType;
    string line;
    string par;

    parFileName = fileName;

    cout << "Reading input from file " << parFileName << endl;

    strcpy(inputfile, parFileName.c_str());
    ifstream myfile(inputfile);

    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;
	if (par == "ParType")
	    {
	    iss >> parType;
	    }
	}
    if (parType == "Positional")
	{
	fileType = parType;
	return 0;
	}
    else if (parType == "Orbital")
	{
	fileType = parType;
	return 1;
	}
    else
	{
	return 2;
	}
    }

int parFile::readParFile()
    {

    int type = parType();

    if (type == 0)
	{
	printf("Read Positions \n");
	readPosFile();

	}
    else if (type == 1)
	{
	printf("Read orbital parameters \n");
	readOrbFile();

	}
    else if (type > 1)
	{
	cout << "Missing header or no parameter file" << endl;

	}
    return type;
    }

int parFile::readParFile(string filename)
    {

    int type = parType(filename);

    if (type == 0)
	{
	printf("Read Positions \n");
	readPosFile();
	}
    else if (type == 1)
	{
	printf("Read orbital parameters \n");
	readOrbFile();

	}
    else if (type > 1)
	{
	cout << "Missing header or no parameter file" << endl;
	}
    return type;
    }

