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

double solradToAU = 0.00467580213904;

parFile::parFile()
    {
    totalMass = 0.0;

    restart = false;
    illumination = false;
    tidal = false;

    systemTime = 0.0;
    maximumTime = 0.0;
    snapshotTime = 0.0;
    snapshotNumber =0;
    nPoints = 0;
    number_bodies = 0;

    }

parFile::parFile(string name)
    {
    parFileName = name;
    totalMass = 0.0;

    restart = false;
    illumination = false;
    tidal = false;

    systemTime = 0.0;
    maximumTime = 0.0;
    snapshotTime = 0.0;
    snapshotNumber = 0;

    nPoints =0;
    number_bodies = 0;

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
    string BodyType,BodyName, meltChoice, restartChoice, illumChoice, tidalChoice;

    int bodyIndex;

    char inputfile[100];

    double val_i, val_j, val_k;

    NBodyFile = "nbody_output.txt";
    snapshotNumber = 0;
    nPoints = 0;

    restart = false;
    tidal = false;
    illumination = false;

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

	if (par == "Restart")
	    {
	    iss >> restartChoice;
	    if (restartChoice == "T")
		{
		restart = true;
		cout << "This run is a restart" << endl;
		}
	    }

	if (par == "PlanetaryIllumination")
	    {
	    iss >> illumChoice;
	    if (illumChoice == "T")
		{
		illumination = true;
		cout << "Planetary Illumination Active" << endl;
		}
	    }

	if (par == "TidalHeating")
	    {
	    iss >> tidalChoice;
	    if (tidalChoice == "T")
		{
		tidal = true;
		cout << "Tidal Heating Active" << endl;
		}
	    }

	if (par == "NBodyOutput" or par == "OutputFile")
	    {
	    iss >> NBodyFile;
	    }
	if (par == "SnapshotTime")
	    {
	    iss >> snapshotTime;
	    }

	if (par == "MaximumTime")
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

	    // Set up vectors to hold all data

	    BodyNames.assign(number_bodies,"");
	    BodyTypes.assign(number_bodies,"");

	    Mass.assign(number_bodies, 0.0);
	    Radius.assign(number_bodies, 0.0);

	    x_position.assign(number_bodies,0.0);
	    y_position.assign(number_bodies,0.0);
	    z_position.assign(number_bodies,0.0);

	    x_velocity.assign(number_bodies,0.0);
	    y_velocity.assign(number_bodies,0.0);
	    z_velocity.assign(number_bodies,0.0);

	    luminosity.assign(number_bodies, 0.0);
	    albedo.assign(number_bodies,0.0);


	    rotationPeriod.assign(number_bodies, 0.0);
	    obliquity.assign(number_bodies, 0.0);
	    winterSolstice.assign(number_bodies, 0.0);
	    oceanFraction.assign(number_bodies, 0.0);
	    initialTemperature.assign(number_bodies, 0.0);
	    orbitCentre.assign(number_bodies, 0.0);
	    activateMelt.assign(number_bodies,false);
	    bodyIndex = -1;

	    }
	else if(par =="NGridPoints")
	    {
	    iss >> nPoints;
	    }
	else if(par=="BodyName")
	    {
	    iss >> BodyName;
	    bodyIndex++;
	    BodyNames[bodyIndex] = BodyName;

	    }
	else if (par == "BodyType")
	    {
	    iss >> BodyType;
	    BodyTypes[bodyIndex]=BodyType;
	    }
	else if (par == "Mass")
	    {
	    iss >> val_i;
	    Mass[bodyIndex]=val_i;
	    totalMass += val_i;
	    }
	else if (par == "Radius")
	    {
	    iss >> val_i;
	    Radius[bodyIndex]=val_i*solradToAU;
	    }
	else if (par == "Position")
	    {
	    iss >> val_i >> val_j >> val_k;
	    x_position[bodyIndex]=val_i;
	    y_position[bodyIndex]=val_j;
	    z_position[bodyIndex]=val_k;
	    }
	else if (par == "Velocity")
	    {
	    iss >> val_i >> val_j >> val_k;
	    x_velocity[bodyIndex]=val_i;
	    y_velocity[bodyIndex]=val_j;
	    z_velocity[bodyIndex]=val_k;
	    }

	else if (par == "Luminosity")
	    {
	    iss >> val_i;
	    luminosity[bodyIndex]=val_i;

	    }

	else if (par == "Albedo")
	    {
	    iss >> val_i;
	    albedo[bodyIndex]=val_i;

	    }

	else if (par == "RotationPeriod")
	    {
	    iss >> val_i;
	    rotationPeriod[bodyIndex]=val_i;

	    }

	else if (par == "Obliquity")
	    {
	    iss >> val_i;
	    val_i = val_i *3.1415/180.0;
	    obliquity[bodyIndex]=val_i;

	    }

	else if (par == "WinterSolstice")
	    {
	    iss >> val_i;
	    val_i = val_i *3.1415/180.0;
	    winterSolstice[bodyIndex]=val_i;

	    }

	else if (par == "OceanFraction")
	    {
	    iss >> val_i;
	    oceanFraction[bodyIndex] = val_i;

	    }
	else if (par == "InitialTemperature")
	    {
	    iss >> val_i;
	    initialTemperature[bodyIndex] = val_i;
	    }
	else if(par == "IceMeltingOn")
	    {
	    iss >> meltChoice;
	    if(meltChoice=="T")
		{
		activateMelt[bodyIndex] = true;
		}

	    }

	}
    myfile.close();

    systemTime = 0.0;
    if(restart)
	{
	setupRestartPositions();
	}

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

    int bodyIndex;
    string par;
    string line;
    string BodyType, BodyName, meltChoice,restartChoice, illumChoice, tidalChoice;

    NBodyFile = "nbody_output.txt";
    snapshotNumber = 0;
    restart = false;
    tidal = false;
    illumination = false;

    char inputfile[100];

    double val_i;

    strcpy(inputfile, parFileName.c_str());

    ifstream myfile(inputfile);


    nPoints = 0;

    // Then loop through each line using getline and then
    //assign to vectors
    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;

	if(par=="Restart")
	    {
	    iss >> restartChoice;
	    if(restartChoice == "T")
		{
		restart=true;
		cout << "This run is a restart" << endl;
		}
	    }

	if (par == "PlanetaryIllumination")
	    {
	    iss >> illumChoice;
	    if (illumChoice == "T")
		{
		illumination = true;
		cout << "Planetary Illumination Active" << endl;
		}
	    }

	if (par == "TidalHeating")
	    {
	    iss >> tidalChoice;
	    if (tidalChoice == "T")
		{
		tidal = true;
		cout << "Tidal Heating Active" << endl;
		}
	    }

	if(par=="NBodyOutput")
	    {
	    iss >> NBodyFile;
	    }
	if(par=="SnapshotTime")
	    {
	    iss >> snapshotTime;
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

	    // Set up vectors to hold all data

	    BodyNames.assign(number_bodies,"");
	    BodyTypes.assign(number_bodies,"");
	    Mass.assign(number_bodies, 0.0);
	    Radius.assign(number_bodies, 0.0);
	    semiMajorAxis.assign(number_bodies, 0.0);
	    eccentricity.assign(number_bodies, 0.0);
	    inclination.assign(number_bodies, 0.0);

	    longAscend.assign(number_bodies, 0.0);
	    inclination.assign(number_bodies, 0.0);
	    Periapsis.assign(number_bodies, 0.0);
	    meanAnomaly.assign(number_bodies, 0.0);

	    luminosity.assign(number_bodies, 0.0);
	    albedo.assign(number_bodies,0.0);

	    rotationPeriod.assign(number_bodies, 0.0);
	    obliquity.assign(number_bodies, 0.0);
	    winterSolstice.assign(number_bodies, 0.0);
	    oceanFraction.assign(number_bodies, 0.0);
	    initialTemperature.assign(number_bodies, 0.0);
	    orbitCentre.assign(number_bodies,0.0);
	    activateMelt.assign(number_bodies,false);
	    bodyIndex = -1;


	    }

	else if(par =="NGridPoints")
		    {
		    iss >> nPoints;
		    }

	else if (par == "BodyName")
	    {
	    iss >> BodyName;
	    bodyIndex++;
	    BodyNames[bodyIndex] = BodyName;
	    }

	else if (par == "BodyType")
	    {
	    iss >> BodyType;
	    BodyTypes[bodyIndex] = BodyType;
	    }
	else if (par == "Mass")
	    {
	    iss >> val_i;
	    Mass[bodyIndex] = val_i;
	    totalMass += val_i;

	    }
	else if (par == "Radius")
	    {
	    // Radii given in solar radii --> convert them to AU
	    iss >> val_i;
	    Radius[bodyIndex] = val_i*solradToAU;
	    }
	else if (par == "SemiMajorAxis")
	    {
	    iss >> val_i;
	    semiMajorAxis[bodyIndex] = val_i;
	    }
	else if (par == "Eccentricity")
	    {
	    iss >> val_i;
	    eccentricity[bodyIndex] = val_i;
	    }
	else if (par == "Inclination")
	    {
	    iss >> val_i;
	    inclination[bodyIndex] = val_i;
	    }
	else if (par == "LongAscend")
	    {
	    iss >> val_i;
	    longAscend[bodyIndex] = val_i;
	    }

	else if (par == "Periapsis")
	    {
	    iss >> val_i;
	    Periapsis[bodyIndex] = val_i;

	    }
	else if (par == "MeanAnomaly")
	    {
	    iss >> val_i;
	    meanAnomaly[bodyIndex] = val_i;

	    }
	else if(par == "OrbitCentre")
	    {

	    iss >> val_i;
	    orbitCentre[bodyIndex] = val_i;
	    }

	else if (par == "Luminosity")
	    {
	    iss >> val_i;
	    luminosity[bodyIndex] = val_i;
	    }

	else if (par == "Albedo")
	    {
	    iss >> val_i;
	    albedo[bodyIndex]=val_i;
	    }

	else if (par == "RotationPeriod")
	    {
	    iss >> val_i;
	    rotationPeriod[bodyIndex] = val_i;

	    }

	else if (par == "Obliquity")
	    {
	    iss >> val_i;
	    val_i = val_i *3.1415/180.0;
	    obliquity[bodyIndex] = val_i;
	    }

	else if (par == "WinterSolstice")
	    {
	    iss >> val_i;
	    val_i = val_i *3.1415/180.0;
	    winterSolstice[bodyIndex] = val_i;
	    }

	else if (par == "OceanFraction")
	    {
	    iss >> val_i;

	    oceanFraction[bodyIndex] = val_i;
	    }
	else if (par == "InitialTemperature")
	    {
	    iss >> val_i;
	    initialTemperature[bodyIndex] = val_i;
	    }
	else if(par == "IceMeltingOn")
	    {
	    iss >> meltChoice;
	    if(meltChoice=="T")
		{
		activateMelt[bodyIndex] = true;
		}

	    }
	}
    myfile.close();

    systemTime = 0.0;

    cout << BodyNames[0] << endl;
    if(restart)
	{
	setupRestartPositions();
	}


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

void parFile::setupRestartPositions()
    {
    /*
     * Written 13/5/14 by dh4gan
     * If the system is restarting from a previous dump, this method
     * generates the correct starting positions and velocities
     */

    int numLines = 0;
    int ibody;
    double blank;
    string line, name;

    // Read Last N Lines of NBody File for time, orbital elements

    cout << "Generating Positions for system restart " << endl;
    cout << "Reading input from file " << NBodyFile << endl;
    cout << "Attempting data read for " << number_bodies << " bodies " << endl;

    ifstream myfile(NBodyFile.c_str());

    if (myfile)
	{

	// First, get number of lines in total

	while (getline(myfile, line))
	    {
	    numLines++;
	    }

	cout << "N Body output file has " << numLines << " lines " << endl;

	}
    else
	{
	cout << "Warning: File " << NBodyFile << " not found" << endl;
	cout << "Assuming that is not a restart!" << endl;
	restart=false;
	return;
	}

    myfile.close();

    // Now read final N lines


    Mass.assign(number_bodies, 0.0);
    Radius.assign(number_bodies, 0.0);

    x_position.assign(number_bodies, 0.0);
    y_position.assign(number_bodies, 0.0);
    z_position.assign(number_bodies, 0.0);

    x_velocity.assign(number_bodies, 0.0);
    y_velocity.assign(number_bodies, 0.0);
    z_velocity.assign(number_bodies, 0.0);

    myfile.open(NBodyFile.c_str());

    int iline = 0;
    while (getline(myfile, line))
	{
	iline++;

	// If at last lines of file, then read information
	if (iline > numLines - number_bodies)
	    {

	    ibody = numLines - iline;
	    ibody = number_bodies - ibody - 1;

	    // Strip commas from line
	    while (line.find(",") != line.npos)
		{
		line.replace(line.find(","), 1, " ");

		}
	    istringstream iss(line);

	    iss >> systemTime;
	    iss >> blank;
	    iss >> name;

	    if (name != BodyNames[ibody])
		{
		cout << "WARNING! Body Names Mismatch: " << name << "  "
			<< BodyNames[ibody] << endl;
		}

	    // Mass, Radius
	    iss >> Mass[ibody];
	    iss >> Radius[ibody];

	    // X, Y, Z, VX, VY, VZ

	    iss >> x_position[ibody];
	    iss >> y_position[ibody];
	    iss >> z_position[ibody];

	    iss >> x_velocity[ibody];
	    iss >> y_velocity[ibody];
	    iss >> z_velocity[ibody];

	    // Orbital Parameters

	    iss >> semiMajorAxis[ibody];
	    iss >> eccentricity[ibody];
	    iss >> inclination[ibody];
	    iss >> longAscend[ibody];
	    iss >> Periapsis[ibody];
	    iss >> meanAnomaly[ibody];

	    }
	}

    }
