/*
 * parFile.cpp
 *
 *  Created on: Sep 23, 2013
 *
 *      Author: davidharvey
 */

#include "parFile.h"
#include "Constants.h"
#include "CSCycleConstants.h"
#include <iostream> // Included for debug lines only
#include <math.h>
#include <string>


parFile::parFile()
    {
   setVariableLocations();

    }

parFile::parFile(string name)
    {
    parFileName = name;
    
    setVariableLocations();
        
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
     * Extracts a Vector3D object containing the position of body index
     *
     */

    Vector3D velocity;

    velocity.elements[0] = vectorDoubleVariables["XVelocity"][index];
    velocity.elements[1] = vectorDoubleVariables["YVelocity"][index];
    velocity.elements[2] = vectorDoubleVariables["ZVelocity"][index];

    return velocity;

    }


void parFile::setVariableType(const vector<const string> &variables, const string &type)

{
    
    int nEntries = variables.size();
    
    for (int i=0;i<nEntries; i++)
    {
        variableLocations[variables[i]] = type;
        //printf("%s %s \n",variables[i].c_str(), type.c_str());
    }
    
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


void parFile::readVariable(string &par, istringstream &iss, int &bodyIndex)

{
    string value;
    //iss >> value;
    
   
    printf("Reading %s \n", par.c_str());

    if(variableLocations[par]==stringType) {readStringVariable(par,iss);}
    
    else if(variableLocations[par]==intType){
        readIntVariable(par,iss);
    
        // If we have read Number of Bodies, then initialise vectors
        if(par=="Number_Bodies") {
            printf("Number of Bodies read as %i \n", intVariables["Number_Bodies"]);
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
            printf("BodyIndex: %i\n",bodyIndex);
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
            printf("Total mass at %i, %f",bodyIndex, doubleVariables["TotalMass"]);
        }
        
    }
    
    else
    {
        cout << "ERROR: parameter " << par << " not recognised " << endl;
    }
    
    
}


void parFile::read3DVector(string &par,istringstream &iss,int &bodyIndex)

{
    double x,y,z;
    iss >> x >> y >> z;

    printf("3D VECTOR %f %f %f \n", x,y,z);
    doubleVariables["X"+par] = x;
    doubleVariables["X"+par] = y;
    doubleVariables["X"+par] = z;
    
}

void parFile::readStringVariable(string &par,istringstream &iss)

{
    string value;
    iss >> value;
    stringVariables[par] = value;
}

void parFile::readIntVariable(string &par,istringstream &iss)

{
    int value;
    iss >> value;
    
    printf("Reading: %s, %i \n", par.c_str(),value);
    intVariables[par] = value;
}

void parFile::readDoubleVariable(string &par,istringstream &iss)

{
    double value;
    iss >> value;
    doubleVariables[par] = value;
}

void parFile::readVectorIntVariable(string &par, istringstream &iss, int &bodyIndex)
{
    int value;
    iss>> value;
    
    vectorIntVariables[par][bodyIndex] = value;
}

void parFile::readVectorDoubleVariable(string &par, istringstream &iss, int &bodyIndex)
{
    double value;
    iss>> value;
    
    vectorDoubleVariables[par][bodyIndex] = value;
}

void parFile::readVectorStringVariable(string &par, istringstream &iss, int &bodyIndex)
{
    string value;
    iss>> value;
    
    vectorStringVariables[par][bodyIndex] = value;
    
    printf("Read vector string variable %s \n", vectorStringVariables[par][bodyIndex].c_str());
}


void parFile::initialiseVectors(int nBodies)
{
    /*
     * Written 19/09/2018 by dh4gan
     * Once the number of bodies is known, creates storage to hold all data
     */
    
    // Assign zeros to int vectors
    
    printf("Initialising vectors %i\n",nBodies);
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

    
    printf("BodyNames size: %i \n", int(vectorStringVariables["BodyName"].size()));
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
    
    
    // If Carbonate Silicate Cycle is on, are stellar spectral types read in?
    // Are Weathering coefficents non-zero?
    
}


void parFile::displayParameters()

{
    
    printf("Global Parameters: \n");
    printf("*********************\n");
    printf("System Name: %s \n",stringVariables["SystemName"].c_str());
    printf("Number of Bodies: %i \n",intVariables["Number_Bodies"]);
    printf("N Body Output File: %s \n",stringVariables["NBodyOutput"].c_str());
    printf("Maximum Time: %f years \n", doubleVariables["MaximumTime"]);
    printf("Snapshot Time: %f years \n", doubleVariables["SnapshotTime"]);
    
    
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
    
    
    printf("Individual Body Parameters \n");
    printf("*********************\n");
    
    for (int i=0; i<intVariables["Number_Bodies"]; i++)
        
    {
        
        if(stringVariables["ParType"].compare("Positional")==0)
        {
        // Write Type, Position and Velocity
        printf("Body %i: Name %s, Type %s \n",i,vectorStringVariables["BodyName"][i].c_str(),vectorStringVariables["BodyType"][i].c_str());
        printf("Position: \n");
        getBodyPosition(i).printVector();
        
        printf("Velocity: \n");
        getBodyVelocity(i).printVector();
    }
        
        else if(stringVariables["ParType"].compare("Orbital")==0)
        {
            printf("Body %i: Name %s, Type %s \n",i,vectorStringVariables["BodyName"][i].c_str(),vectorStringVariables["BodyType"][i].c_str());
            printf("Orbit: a e i LongAscend MeanAnomaly\n");
            printf("%f %f %f %f %f %f \n",vectorDoubleVariables["SemiMajorAxis"][i], vectorDoubleVariables["Eccentricity"][i],vectorDoubleVariables["Inclination"][i],vectorDoubleVariables["LongAscend"][i],vectorDoubleVariables["Periapsis"][i],vectorDoubleVariables["MeanAnomaly"][i]);
        }
        }
    
    
}

void parFile::convertToRadians(int nBodies)

{
    
    string degreeVariables[] = {"Obliquity","WinterSolstice"};
    
    for (int i=0; i< sizeof(degreeVariables)/sizeof(*degreeVariables); i++)
    {
        for (int ibody=0; ibody<nBodies; ibody++)
        {
        vectorDoubleVariables[degreeVariables[i]][ibody] = vectorDoubleVariables[degreeVariables[i]][ibody]*degToRad;
        
    }
    }
    
}

bool parFile::initialiseBoolean(string &par)

{
    bool choice=false;
    if(stringVariables[par].compare("T")==0 or stringVariables[par].compare("t")==0 or stringVariables[par].compare("Y")==0 or stringVariables[par].compare("y")==0 )
    {
        boolVariables[par]=true;
    }
    
    return choice;
}


void parFile::initialiseAllBooleans()

{
    // Booleans
    string boolVar[] = {"Restart", "ObliquityEvolution","CarbonateSilicateCycle","TidalHeating","PlanetaryIllumination","FullOutput"};
    
    for (int i=0; i<sizeof(boolVar); i++)
    {
        boolVariables[boolVar[i]] = initialiseBoolean(boolVar[i]);
    }
    
}



void parFile::readFile()

{

    cout << "What is the input file? " << endl;

    getline(cin, parFileName);

    cout << "Reading input from file " << parFileName << endl;
    
    parFile(filename);
}




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
        istringstream iss(line);
        iss >> par;
        
        if(par.compare("--")!=0)
        {
        readVariable(par,iss,bodyIndex);
        }
        
    }
    myfile.close();
    
    printf("File Read complete\n");
    
    // Convert any variables read in degrees to radians
    convertToRadians(intVariables["Number_Bodies"]);
    printf("Quantities converted to Radians\n");
    
    doubleVariables["SystemTime"] = 0.0;
    if(boolVariables["Restart"])
    {
        setupRestartPositions();
    }
    printf("System time is %f\n", doubleVariables["SystemTime"]);
}

/*
void parFile::readPosFile()
    {
    /*
     * Read par files from input and set up the body and parse back
     * body
     *
     */
 /*   // First read in the user data file
    // Format must be,name of body followed by its parameters

    string par;
    string line;
    string BodyType,BodyName,spectral;
    string meltChoice, restartChoice, illumChoice, tidalChoice;
    string obliqChoice, CScycleChoice, fullOutputChoice;

    int bodyIndex=-1;

    char inputfile[100];

    double val_i, val_j, val_k;

    NBodyFile = "nbody_output.txt";
    snapshotNumber = 0;
    nPoints = 0;

    restart = false;

    tidalHeatingOn = false;
    illuminationOn = false;
    obliquityOn = false;

    CSCycleOn = false; //Giblin 10/7/15

    strcpy(inputfile, parFileName.c_str());

    ifstream myfile(inputfile);

    // Then loop through each line using getline and then
    //assign to vectors

    // Handy defaults for incomplete files

    luminosity.push_back(0.0);
    rotationPeriod.push_back(0.0);
    obliquity.push_back(0.0);
    precession.push_back(0.0);
    oceanFraction.push_back(0.0);
    initialTemperature.push_back(0.0);

    while (getline(myfile, line))
	{
	istringstream iss(line);
	iss >> par;
        
        //readVariable(par,iss,bodyIndex);

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
		illuminationOn = true;
		cout << "Planetary Illumination Active" << endl;
		}
	    }

	if (par == "TidalHeating")
	    {
	    iss >> tidalChoice;
	    if (tidalChoice == "T")
		{
		tidalHeatingOn = true;
		cout << "Tidal Heating Active" << endl;
		}
	    }

	if (par == "ObliquityEvolution")
	    {
	    iss >> obliqChoice;
	    if (obliqChoice == "T")
		{
		obliquityOn = true;
		cout << "Obliquity Evolution Active" << endl;
		}
	    }


	if (par == "CarbonateSilicateCycle")
	    {
	    iss >> CScycleChoice;
	    if (CScycleChoice == "T")
		{
		CSCycleOn = true;
		cout << "Carbonate Silicate Cycle Active" << endl;
		}
	    }


	if(par == "FullOutput")
	{
		iss >> fullOutputChoice;

		if(fullOutputChoice =="y")
		{
			fullOutput=true;
			cout<<"Full output ON LADDDDDDDDD"<<endl; //debugging
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
	    spectralType.assign(number_bodies, "G"); // G is the default spectral type
	    albedo.assign(number_bodies,0.0);


	    rotationPeriod.assign(number_bodies, 0.0);
	    obliquity.assign(number_bodies, 0.0);
	    precession.assign(number_bodies, 0.0);
	    oceanFraction.assign(number_bodies, 0.0);
	    initialTemperature.assign(number_bodies, 0.0);
	    outgassingRate.assign(number_bodies,outgassingRateEarth); // outgassing rate set to Earth by default
	    seafloorWeathering.assign(number_bodies,0.0); // seafloor weathering off by default
	    betaCO2.assign(number_bodies, 0.5); // Set this to 0.5 as default
	    gammaCO2.assign(number_bodies,1.0); // Set this to unity as default
	    orbitCentre = vector<int>(number_bodies,0);
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

	else if (par == "SpectralType")
	    {

	    iss >> spectral;
	    spectralType[bodyIndex] = spectral;

	    }

	else if (par == "Albedo") // For Planet objects only
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
	    val_i = val_i *degToRad;
	    obliquity[bodyIndex]=val_i;

	    }

	else if (par == "Precession" or par =="WinterSolstice")
	    {
	    iss >> val_i;
	    val_i = val_i *degToRad; // Convert to radians
	    precession[bodyIndex]=val_i;

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
	else if (par == "outgassingRate")
	    {
	    iss >> val_i;
	    outgassingRate[bodyIndex] = val_i;
	    }
	else if (par == "landWeatheringParameter")
	    {
	    iss >> val_i;
	    betaCO2[bodyIndex] = val_i;
	    }
	else if(par == "oceanWeatheringRate")
	    {
	    iss >> val_i;
	    seafloorWeathering[bodyIndex] = val_i;
	    }

	else if(par =="abioticWeatheringParameter")
	    {
	    iss >> val_i;
	    betaCO2[bodyIndex] = val_i;
	    }
	else if(par =="oceanWeatheringParameter")
	    {
	    iss >> val_i;
	    gammaCO2[bodyIndex] = val_i;
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

/*    int bodyIndex;
    string par;
    string line;
    string BodyType, BodyName,spectral;

    string meltChoice, restartChoice, illumChoice, tidalChoice;
    string obliqChoice, CScycleChoice, fullOutputChoice;

    NBodyFile = "nbody_output.txt";
    snapshotNumber = 0;
    restart = false;
    tidalHeatingOn = false;
    illuminationOn = false;
    obliquityOn = false;

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
		illuminationOn = true;
		cout << "Planetary Illumination Active" << endl;
		}
	    }

	if (par == "TidalHeating")
	    {
	    iss >> tidalChoice;
	    if (tidalChoice == "T")
		{
		tidalHeatingOn = true;
		cout << "Tidal Heating Active" << endl;
		}
	    }

	if (par == "ObliquityEvolution")
	    {
	    iss >> obliqChoice;
	    if (obliqChoice == "T")
		{
		obliquityOn = true;
		cout << "Obliquity Evolution Active" << endl;
		}
	    }


	if (par == "CarbonateSilicateCycle")
	    {
	    iss >> CScycleChoice;
	    if (CScycleChoice == "T")
		{
		CSCycleOn = true;
		cout << "Carbonate Silicate Cycle Active" << endl;
		}
	    }

	if(par == "FullOutput")
	{
		iss >> fullOutputChoice;

		if(fullOutputChoice =="y")
		{
			fullOutput=true;
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
	    spectralType.assign(number_bodies, "G"); // G is the default spectral type
	    albedo.assign(number_bodies,0.0);

	    rotationPeriod.assign(number_bodies, 0.0);
	    obliquity.assign(number_bodies, 0.0);
	    precession.assign(number_bodies, 0.0);
	    oceanFraction.assign(number_bodies, 0.0);
	    initialTemperature.assign(number_bodies, 0.0);
	    outgassingRate.assign(number_bodies,outgassingRateEarth);
	    seafloorWeathering.assign(number_bodies,0.0); // seafloor weathering off by default
	    betaCO2.assign(number_bodies, 0.5); // Set this to 0.5 as default
	    gammaCO2.assign(number_bodies,1.0); // Set this to unity as default


	    orbitCentre = vector<int>(number_bodies,0);
	    activateMelt.assign(number_bodies,false);
	    bodyIndex = -1;

	    //possibly require a line here for CO2Pressure vector...?

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
	    orbitCentre[bodyIndex] = int(val_i);
	    }

	else if (par == "Luminosity")
	    {
	    iss >> val_i;
	    luminosity[bodyIndex] = val_i;
	    }

	else if (par == "SpectralType")
	    {

	    iss >> spectral;
	    spectralType[bodyIndex] = spectral;

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
	    val_i = val_i *degToRad;
	    obliquity[bodyIndex] = val_i;
	    }

	else if (par == "Precession" or par =="WinterSolstice")
	    {
	    iss >> val_i;
	    val_i = val_i *degToRad;
	    precession[bodyIndex] = val_i;
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
	else if (par == "outgassingRate")
	    {
	    iss >> val_i;
	    outgassingRate[bodyIndex] = val_i;
	    }
	else if (par == "landWeatheringParameter")
	    {
	    iss >> val_i;
	    betaCO2[bodyIndex] = val_i;
	    }
	else if(par == "oceanWeatheringRate")
	    {
	    iss >> val_i;
	    seafloorWeathering[bodyIndex] = val_i;
	    }
	else if(par =="oceanWeatheringParameter")
	    {
	    iss >> val_i;
	    gammaCO2[bodyIndex] = val_i;
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


    }*.
/*int parFile::parType()
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
/*    char inputfile[100];
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
 /*   char inputfile[100];
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
*/
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

	    iss >> doubleVariables["systemTime"];
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

	    iss >> vectorDoubleVariables["semiMajorAxis"][ibody];
	    iss >> vectorDoubleVariables["eccentricity"][ibody];
	    iss >> vectorDoubleVariables["Inclination"][ibody];
	    iss >> vectorDoubleVariables["LongAscend"][ibody];
	    iss >> vectorDoubleVariables["Periapsis"][ibody];
	    iss >> vectorDoubleVariables["MeanAnomaly"][ibody];

	    }
	}

    }
