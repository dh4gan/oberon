/*
 * parFile.cpp
 *
 *  Created on: Sep 23, 2013
 *  Majorly revised: Sep 19, 2018
 *
 *      Author: dh4gan
 */

#include "parFile.h"
#include "Constants.h"
#include "CSCycleConstants.h"
#include <iostream> // Included for debug lines only
#include <math.h>
#include <string>


// Constructors for the parFile object

parFile::parFile()
    {
   setVariableLocations();

    }

parFile::parFile(string name)
    {
    parFileName = name;
    
    setVariableLocations();
        
    }

parFile::~parFile() {
}

// End of constructors

void parFile::readFile(string &filename)

{
    /*
     * Written 20/09/18 by dh4gan
     * Read in data from the OBERON parameter file
     *
     */
    
    string par;
    string line;
    
    int bodyIndex=-1;
    
    ifstream myfile(filename.c_str());
    
    // Then loop through each line using getline and then
    //assign to vectors
    
    while (getline(myfile, line))
    {
	par ="";
        istringstream iss(line);
        iss >> par;

	if(par.compare("")==0) break;        
        if(par.compare("--")!=0)
        {
            readVariable(par,iss,bodyIndex);
        }
        
    }
    myfile.close();
    
    printf("\t\tFile Read complete\n");
    printf("%s",screenBar.c_str());
    
    // Convert any variables read in degrees to radians
    convertToRadians(intVariables["Number_Bodies"]);
    
    // Set up boolean variables
    initialiseAllBooleans();
    
    doubleVariables["SystemTime"] = 0.0;
    
    if(boolVariables["Restart"])
    {
        setupRestartPositions();
    }
}


void parFile::readFile()

{
    /*
     *
     * Written 20/09/2019 by dh4gan
     * Overloaded method for readFile(filename)
     * allows the user to enter the parameter filename at the command line
     * before reading in
     *
     */
    
    cout << "What is the input file? " << endl;
    
    getline(cin, parFileName);
    
    cout << "Reading input from file " << parFileName << endl;
    
    parFile(filename);
    
}


void parFile::readVariable(string &par, istringstream &iss, int &bodyIndex)

{
    /*
     * Written 20/09/2018 by dh4gan
     * Given the input parameter in question, it reads data from the parameter file into the appropriate map
     * Inputs:
     *
     * par - string describing the parameter
     * iss - stringstream object doing the file read
     * bodyIndex - which body are we currently reading data for?
     *
     */
    
    
    string value;
    
    if(variableLocations[par]==stringType) {readStringVariable(par,iss);}
    
    else if(variableLocations[par]==intType){
        readIntVariable(par,iss);
        
        // If we have read Number of Bodies, then initialise vectors
        if(par=="Number_Bodies") {
            initialiseVectors(intVariables["Number_Bodies"]);
        }
    }
    
    else if(variableLocations[par]==doubleType){
        
        readDoubleVariable(par,iss);
        
    }
    else if(variableLocations[par]==vectorStringType){
        
        // If reading BodyName, assume we have moved on to the next body
        if(par=="BodyName")
        {
            bodyIndex++;
        }
        
        readVectorStringVariable(par,iss,bodyIndex);
    }
    else if(variableLocations[par]==vectorIntType){readVectorIntVariable(par,iss,bodyIndex);}
    else if(variableLocations[par]==vectorDoubleType)
    {
        
        // If reading 3D vectors, ensure these are stored correctly
        if(par=="Position" or par=="Velocity")
        {
            read3DVector(par,iss,bodyIndex);
        }
        else
        {
            readVectorDoubleVariable(par,iss,bodyIndex);
        }
        
        if(par.compare("Mass")==0)
        {
            doubleVariables["TotalMass"] = doubleVariables["TotalMass"] + vectorDoubleVariables["Mass"][bodyIndex];
        }
        
        if(par.compare("OutgassingRate")==0)
        {
            vectorDoubleVariables["OutgassingRate"][bodyIndex] *= outgassingRateEarth;
        }
        
    }
    
    else
    {
        cout << "ERROR: parameter " << par << " not recognised " << endl;
    }
    
    
}



Vector3D parFile::getBodyPosition(int index)
    {
    /*
     * Author: dh4gan
     * Extracts a Vector3D object containing the position of body index
     *
     */

    Vector3D position;

    position.elements[0] = vectorDoubleVariables["XPosition"][index];
    position.elements[1] = vectorDoubleVariables["YPosition"][index];
    position.elements[2] = vectorDoubleVariables["ZPosition"][index];

    return position;

    }

Vector3D parFile::getBodyVelocity(int index)
    {
    /*
     * Author: dh4gan
     * Extracts a Vector3D object containing the velocity of body index
     *
     */

    Vector3D velocity;

    velocity.elements[0] = vectorDoubleVariables["XVelocity"][index];
    velocity.elements[1] = vectorDoubleVariables["YVelocity"][index];
    velocity.elements[2] = vectorDoubleVariables["ZVelocity"][index];

    return velocity;

    }



    

void parFile::setVariableLocations()

{
    
    /*
     * Written 19/09/2018 by dh4gan
     * Sets up the mapping of variables to map objects
     *
     */
    
    setVariableType(stringVar, stringType);             // string
    setVariableType(boolVar, stringType);                 // bool
    setVariableType(intVar, intType);                   // int
    setVariableType(doubleVar, doubleType);             // double
    setVariableType(vectorStringVar, vectorStringType); // vector<string>
    setVariableType(vectorIntVar, vectorIntType);       // vector<int>
    setVariableType(vectorDoubleVar, vectorDoubleType); // vector<double>

    
}

void parFile::setVariableType(const vector<string> &variables, const string &type)

{
    /*
     * Written 20/09/2018 by dh4gan
     * Assigns a variable type to a vector array (used by setVariableLocations())
     * This allows us to identify which map a given variable is stored in
     *
     */
    
    
    int nEntries = variables.size();
    
    for (int i=0;i<nEntries; i++)
    {
        variableLocations[variables[i]] = type;
        //printf("%s %s \n",variables[i].c_str(), type.c_str());
    }
    
}



void parFile::read3DVector(string &par,istringstream &iss,int &bodyIndex)

{
    
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Given 3D data in the parameter file (position or velocity),
     * correctly reads it in for decanting into a Vector3D object later
     *
     */
     
     
    double x,y,z;
    iss >> x >> y >> z;

    vectorDoubleVariables["X"+par][bodyIndex] = x;
    vectorDoubleVariables["Y"+par][bodyIndex] = y;
    vectorDoubleVariables["Z"+par][bodyIndex] = z;
    
    
}

void parFile::readStringVariable(string &par,istringstream &iss)

{
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Facilitates reading of string data into map object
     */
    
    string value;
    iss >> value;
    stringVariables[par] = value;
}

void parFile::readIntVariable(string &par,istringstream &iss)

{
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Facilitates reading of int data into map object
     */
    
    int value;
    iss >> value;
    
    intVariables[par] = value;
}

void parFile::readDoubleVariable(string &par,istringstream &iss)

{
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Facilitates reading of double data into map object
     */
    
    double value;
    iss >> value;
    doubleVariables[par] = value;
}

void parFile::readVectorIntVariable(string &par, istringstream &iss, int &bodyIndex)
{
    
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Facilitates reading of int data specific to a body into map object
     */
    
    int value;
    iss>> value;
    
    vectorIntVariables[par][bodyIndex] = value;
}

void parFile::readVectorDoubleVariable(string &par, istringstream &iss, int &bodyIndex)
{
    
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Facilitates reading of double data specific to a body into map object
     */
    
    double value;
    iss>> value;
    
    vectorDoubleVariables[par][bodyIndex] = value;
}

void parFile::readVectorStringVariable(string &par, istringstream &iss, int &bodyIndex)
{
    
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Facilitates reading of string data specific to a body into map object
     */
    
    string value;
    iss>> value;
    
    vectorStringVariables[par][bodyIndex] = value;
    
}


void parFile::initialiseVectors(int nBodies)
{
    /*
     * Written 19/09/2018 by dh4gan
     * Once the number of bodies is known, creates storage to hold all data
     */
    
    // Assign zeros to int vectors
    
    std::map < string, vector<int> >:: iterator iti = vectorIntVariables.begin();
    
    vector<int> emptyIntVector(nBodies,0);
    vector<double> emptyDoubleVector(nBodies,0.0);
    vector<string> emptyStringVector(nBodies,"");
    
    for (int i=0; i<int(vectorIntVar.size()); i++)
    {
        vectorIntVariables[vectorIntVar[i]]=emptyIntVector;
    }
    
    for (int i=0; i<int(vectorDoubleVar.size()); i++)
    {
        vectorDoubleVariables[vectorDoubleVar[i]]=emptyDoubleVector;
    }
    
    for (int i=0; i<int(vectorStringVar.size()); i++)
    {
        vectorStringVariables[vectorStringVar[i]]=emptyStringVector;
    }

    // Set useful default values here
    
    for (int i=0; i<nBodies; i++)
    {
    vectorDoubleVariables["BetaCO2"][i] = betaCO2_0;
    vectorDoubleVariables["GammaCO2"][i] = gammaCO2_0;
    vectorDoubleVariables["OutgassingRate"][i] = outgassingRateEarth;
    vectorDoubleVariables["SeafloorWeatheringRate"][i] = seaFloorWeathering_0;
    }
}


void parFile::reportError(const string &par, double &value)

{
    printf("Error in parameter %s: invalid value %f \n",par.c_str(),value);
    exit(1);
}

void parFile::reportError(const string &par, int &value)

{
    printf("Error in parameter %s: invalid value %i \n",par.c_str(),value);
    exit(1);
}


void parFile::reportError(const string &par, string &value)

{
    printf("Error in parameter %s: invalid value %s \n",par.c_str(),value.c_str());
    exit(1);
}

void parFile::checkParameters()

{
    
    /*
     * Written 20/09/2018 by dh4gan
     *
     * Performs consistency checks on all inputs
     *
     */
    
    // Check number of grid points is defined
    if(intVariables["NGridPoints"]==0)
    {
        reportError("NGridPoints",intVariables["nGridPoints"]);
    }
    
    
    // Set a default name for system if not specified
    if(stringVariables["SystemName"].compare("")==0)
    {
        stringVariables["SystemName"] = "System";
    }
    
    
    // Check snapshot interval is defined
    if(doubleVariables["SnapshotTime"]==0.0)
    {
        reportError("SnapshotTime",doubleVariables["SnapshotTime"]);
    }
    
    // Check maximum time is defined
    if(doubleVariables["MaximumTime"]==0.0)
    {
        reportError("MaximumTime", doubleVariables["MaximumTime"]);
    }
               
}


void parFile::displayParameters()

{
    /*
     * Written 24/09/2018 by dh4gan
     * Writes all inputted parameters to the screen
     *
     */
    
    printf("Global Parameters: \n");
    printf("%s",screenBar.c_str());
    printf("System Name: %s \n",stringVariables["SystemName"].c_str());
    printf("Initial Time: %.1E \n",doubleVariables["SystemTime"]);
    printf("Number of Bodies: %i \n",intVariables["Number_Bodies"]);
    printf("N Body Output File: %s \n",stringVariables["NBodyOutput"].c_str());
    printf("Maximum Time: %.1E years \n", doubleVariables["MaximumTime"]);
    printf("Snapshot Time: %.1E years \n", doubleVariables["SnapshotTime"]);
    
    printf("%s",screenBar.c_str());
    if(boolVariables["Restart"])
    {
        printf("This is a restart - using vector data from pre-existing nbody output file \n");
    }
    
    if(boolVariables["PlanetaryIllumination"])
    {
        printf("Planetary Illumination is ON \n");
    }
    else
    {
        printf("Planetary Illumination is OFF \n");
    }
    
    if(boolVariables["TidalHeating"])
    {
        printf("Tidal Heating is ON \n");
    }
    else
    {
        printf("Tidal Heating is OFF \n");
    }
    
    
    if(boolVariables["CarbonateSilicateCycle"])
    {
        printf("Carbonate Silicate Cycle is ON \n");
    }
    else
    {
        printf("Carbonate Silicate Cycle is OFF \n");
    }
    
    printf("%s",screenBar.c_str());
    printf("Individual Body Parameters \n");
    printf("%s",screenBar.c_str());
    
    for (int i=0; i<intVariables["Number_Bodies"]; i++)
        
    {
        
        printf("Body %i: Name %s, Type %s \n",i,vectorStringVariables["BodyName"][i].c_str(),vectorStringVariables["BodyType"][i].c_str());
        
    
        // Print stellar data in solar masses/radii
        if(vectorStringVariables["BodyType"][i].compare("Star")==0)
        {
        printf("Mass: %.2f solar masses, \nRadius: %f solar radii\n",vectorDoubleVariables["Mass"][i], vectorDoubleVariables["Radius"][i]);
        }
        
        // Otherwise print in earth masses/radii
        else
        {
             printf("Mass: %.2f Earth masses, \nRadius: %f Earth radii\n",vectorDoubleVariables["Mass"][i]*msolToMEarth, vectorDoubleVariables["Radius"][i]*solradToREarth);
        }
        
        
        
        if(stringVariables["ParType"].compare("Positional")==0)
        {
        // Write Type, Position and Velocity
        
        printf("Position: \n");
        getBodyPosition(i).printVector();
        
        printf("Velocity: \n");
        getBodyVelocity(i).printVector();
    }
        
        else if(stringVariables["ParType"].compare("Orbital")==0)
        {
            
            printf("Orbit: a e i LongAscend MeanAnomaly\n");
            printf("%f %f %f %f %f %f \n",vectorDoubleVariables["SemiMajorAxis"][i], vectorDoubleVariables["Eccentricity"][i],vectorDoubleVariables["Inclination"][i],vectorDoubleVariables["LongAscend"][i],vectorDoubleVariables["Periapsis"][i],vectorDoubleVariables["MeanAnomaly"][i]);
            
        }
         printf("%s",screenBar.c_str());
    }
    
    
}

void parFile::convertToRadians(int nBodies)

{
    /*
     * Written 25/09/2018 by dh4gan
     * Input angular variables are all in degrees - must be converted to radians
     *
     */
    
    
    for (int i=0; i<degreeVar.size(); i++)
    {
        for (int ibody=0; ibody<nBodies; ibody++)
        {
        vectorDoubleVariables[degreeVar[i]][ibody] = vectorDoubleVariables[degreeVar[i]][ibody]*degToRad;
        
    }
    }
    
}




void parFile::initialiseAllBooleans()

{
    
    /*
     * Written 25/09/2018 by dh4gan
     * Set up all boolean values
     *
     */
    
    for (int i=0; i<boolVar.size(); i++)
    {
        initialiseBoolean(boolVar[i]);
    }
    
}

void parFile::initialiseBoolean(const string &par)

{
    /*
     * Written 25/09/2018 by dh4gan
     * Set up boolean value for parameter given by "par"
     *
     */
    
    boolVariables[par] = false;
    if(stringVariables[par].compare("T")==0 or stringVariables[par].compare("t")==0 or stringVariables[par].compare("Y")==0 or stringVariables[par].compare("y")==0 )
    {
        boolVariables[par]=true;
    }
    
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
    
        
    int Nbodies = intVariables["Number_Bodies"];
    string NBodyFile = stringVariables["NBodyOutput"];
        
    // Read Last N Lines of NBody File for time, orbital elements

    cout << "Generating Positions for system restart " << endl;
    cout << "Reading input from file " << NBodyFile << endl;
    cout << "Attempting data read for " << Nbodies << " bodies " << endl;

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
	boolVariables["Restart"]=false;
	return;
	}

    myfile.close();

    // Now read final N lines
        
    vectorDoubleVariables["Mass"].assign(Nbodies, 0.0);
    vectorDoubleVariables["Radius"].assign(Nbodies, 0.0);

    vectorDoubleVariables["XPosition"].assign(Nbodies, 0.0);
    vectorDoubleVariables["YPosition"].assign(Nbodies, 0.0);
    vectorDoubleVariables["YPosition"].assign(Nbodies, 0.0);

    vectorDoubleVariables["XVelocity"].assign(Nbodies, 0.0);
    vectorDoubleVariables["YVelocity"].assign(Nbodies, 0.0);
    vectorDoubleVariables["ZVelocity"].assign(Nbodies, 0.0);

    myfile.open(NBodyFile.c_str());

    int iline = 0;
    while (getline(myfile, line))
	{
	iline++;

	// If at last lines of file, then read information
	if (iline > numLines - Nbodies)
	    {

	    ibody = numLines - iline;
	    ibody = Nbodies - ibody - 1;

	    // Strip commas from line
	    while (line.find(",") != line.npos)
		{
		line.replace(line.find(","), 1, " ");

		}
	    istringstream iss(line);

	    iss >> doubleVariables["SystemTime"];
	    iss >> blank;
	    iss >> name;
            

	    if (name.compare(vectorStringVariables["BodyName"][ibody])!=0)
		{
		cout << "WARNING! Body Names Mismatch: " << name << "  "
			<< vectorStringVariables["BodyName"][ibody] << endl;
		}

	    // Mass, Radius
	    iss >> vectorDoubleVariables["Mass"][ibody];
	    iss >> vectorDoubleVariables["Radius"][ibody];

	    // X, Y, Z, VX, VY, VZ

	    iss >> vectorDoubleVariables["XPosition"][ibody];
	    iss >> vectorDoubleVariables["YPosition"][ibody];
	    iss >> vectorDoubleVariables["ZPosition"][ibody];

	    iss >> vectorDoubleVariables["XVelocity"][ibody];
	    iss >> vectorDoubleVariables["YVelocity"][ibody];
	    iss >> vectorDoubleVariables["ZVelocity"][ibody];

	    // Orbital Parameters

	    iss >> vectorDoubleVariables["SemiMajorAxis"][ibody];
	    iss >> vectorDoubleVariables["Eccentricity"][ibody];
	    iss >> vectorDoubleVariables["Inclination"][ibody];
	    iss >> vectorDoubleVariables["LongAscend"][ibody];
	    iss >> vectorDoubleVariables["Periapsis"][ibody];
	    iss >> vectorDoubleVariables["MeanAnomaly"][ibody];

	    }
	}

    }
