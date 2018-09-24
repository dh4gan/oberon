/*
 * main.cpp
 *
 *  Created on: Jan 9, 2014
 *      Author: dh4gan
 *
 *	Reads in parameter files and runs N Body code
 *	where some objects have climate modelling done in tandem
 *	via LEBM modelling
 *
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Constants.h"
#include "System.h"
#include "Star.h"
#include "Planet.h"
#include "World.h"
#include "parFile.h"

#include <fstream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[])
    {

    double G = Gmau;
    double unit2sec = year/twopi;


    int i, fileType;
    int snapshotNumber=0;

    double tStop;
    double timeunit,timeyr;
    double dtsec,dtunit;
    double tSnap, tMax;


    Vector3D body_i_position;
    Vector3D body_i_velocity;

    System nBodySystem;
    vector<Body*> BodyArray;
    vector<Body*> Bodies;

    parFile input;
    FILE * outputfile;

    printf("  \n");
    printf("*********************************************** \n");
    printf("    OBERON - OBliquity and Energy balance Run on Nbody systems \n ");
    printf("    Date Created : 9th January 2014 \n");
    printf("    Current Version: 20th September 2018 \n");
    printf("*********************************************** \n");
    printf("  \n");

    // Read in parameters file

    if (argc == 2)
	{
	string fileString = string(argv[1]);
	fileType = input.readParFile(fileString);
	}
    else
	{
	fileType = input.readParFile();
	if (fileType > 1)
	    {
	    return -1;
	    }
	}

    if(input.getIntVariable("NGridPoints")==0)
	{
	printf("ERROR! LEBM grid points not defined\n Is NGridPoints defined in paramsfile?\n");
	return -1;
	}

    // Record parameter data

    tMax = input.getDoubleVariable("MaximumTime");
    tSnap = input.getDoubleVariable("SnapshotTime");

    if(input.getBoolVariable("Restart") == true)
	{
		cout << "Restart - Using vector data from nbody output" << endl;
	}
    

   if(input.getBoolVariable("CarbonateSilicateCycle"))
	{
		cout << "CS cycle is ON" << endl;
	}

    //First loop through each of the bodies in the system and set them up
    //adding them to the BodyArray
    for (i = 0; i < input.getIntVariable("Number_Bodies"); i++)
	{
	if (fileType == 0 or input.getBoolVariable("Restart"))
	    {

	    cout << "setting up body from vectors" << endl;

	    body_i_position = input.getBodyPosition(i);
	    body_i_velocity = input.getBodyVelocity(i);


	    // If the Body is a Star, add a Star Object

	    if (input.BodyTypes[i] == "Star")
		{
		BodyArray.push_back(
			new Star(input.BodyNames[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.luminosity[i], input.spectralType[i]));
		}

	    // If the Body is a Planet, add a Planet Object
	    if (input.BodyTypes[i] == "Planet")
		{
		BodyArray.push_back(
			new Planet(input.BodyNames[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.albedo[i]));
		}

	    // If the Body is a World, add a World Object and set up LEBM

	    if (input.BodyTypes[i] == "World")
		{
		// Code will halt if initial temperature zero
		// this stops incomplete params files running successfully

		if(input.initialTemperature[i]==0.0)
		    {
		    printf("ERROR in World %s Setup: initial T zero \nIs World fully specified in params file? ", input.BodyNames[i].c_str());
		    return -1;
		    }
		BodyArray.push_back(
			new World(input.BodyNames[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.nPoints, input.obliquity[i],input.rotationPeriod[i], input.precession[i],
				input.oceanFraction[i], input.initialTemperature[i], input.activateMelt[i], input.getBoolVariable("Restart"),
				input.tidalHeatingOn, input.obliquityOn, input.CSCycleOn,input.outgassingRate[i], input.betaCO2[i], input.seafloorWeathering[i], input.gammaCO2[i])); //Fed to parFile.cpp

		if(input.getBoolVariable("Restart"))
		    {
		    cout << "Reading Temperature data for World " << BodyArray.back()->getName() << endl;
		    snapshotNumber = BodyArray.back()->getRestartParameters();

		    if(snapshotNumber==-1)
			{
			printf("ERROR in World %s temperature setup \n",input.BodyNames[i].c_str() );
			return -1;
			}
		    }

		}


	    }
	else if (fileType == 1 and input.getBoolVariable("Restart")==false)
	    {
	    printf("setting up body %s with orbital parameters \n", input.BodyNames[i].c_str());

	    // If the Body is a Star, add a Star Object
	    if (input.BodyTypes[i] == "Star")
		{

		BodyArray.push_back(
			new Star(input.BodyNames[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass, input.luminosity[i], input.spectralType[i]));

		}

	    // If the Body is a Planet, add a Planet Object
	    if (input.BodyTypes[i] == "Planet")
		{
		BodyArray.push_back(
			new Planet(input.BodyNames[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass, input.albedo[i]));

		}

	    // If the Body is a World, add a World Object and set up LEBM
	    if (input.BodyTypes[i] == "World")
		{

		if (input.initialTemperature[i] == 0.0)
		    {
		    printf(
			    "ERROR in World %s Setup: initial T zero \nIs World fully specified in params file? \n",
			    input.BodyNames[i].c_str());
		    return -1;
		    }
		BodyArray.push_back(
			new World(input.BodyNames[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass,input.nPoints, input.obliquity[i],input.rotationPeriod[i], input.precession[i],
				input.oceanFraction[i], input.initialTemperature[i], input.activateMelt[i], input.getBoolVariable("Restart"),
				input.tidalHeatingOn, input.obliquityOn, input.CSCycleOn,
				input.outgassingRate[i], input.betaCO2[i], input.seafloorWeathering[i], input.gammaCO2[i]));



		}

	    }

	printf("body %s set up \n", BodyArray.back()->getName().c_str());

	}

    // Set up System object using BodyArray

    printf("Setting up system %s \n", input.SystemName.c_str());

    nBodySystem = System(input.SystemName, BodyArray);

    // If the System is created from orbital parameters, set up vectors here

    nBodySystem.setHostBodies(input.orbitCentre);

    if(fileType ==1 and input.getBoolVariable("Restart")==false)
	{
	nBodySystem.setupOrbits(input.orbitCentre);
	}

    // Calculate its initial properties
    nBodySystem.calcInitialProperties();

    // Switch Planetary Illumination on/off
    nBodySystem.setIllumination(input.illuminationOn);


    // Set up the outputs

    if (input.getBoolVariable("Restart") and snapshotNumber !=0)
	{
	outputfile = fopen(input.getStringVariable("NBodyOutput"), "a");
	}
    else
	{
	outputfile = fopen(input.getStringVariable("NBodyOutput"), "w");
	fprintf(outputfile, "Number of Bodies, %i \n", input.getIntVariable("Number_Bodies"));
	}


    // Now loop over snap shots, outputting the system data each time

    tStop = 0.0;
    tMax = tMax * twopi; // Convert maximum time to code units
    tSnap = tSnap*twopi; // Convert snapshot time to code units

    // Timesteps will be calculated in NBody units, and converted back to seconds for LEBM

    // Calculate the minimum LEBM timestep for all worlds and NBody timestep
    dtunit = nBodySystem.calcCombinedTimestep();


    // If the system is restarting, ramp down the first timestep for stability
    if(input.getBoolVariable("Restart"))
	{
	dtunit = dtunit/100.0;
	}

    dtsec = dtunit*unit2sec; // This will be the default timestep measure - derive dtunit = dtsec*2pi/(3.15e7)
    timeunit = input.systemTime*twopi;


    printf("System set up: Running \n");
    while (timeunit < tMax)
	{
	tStop = timeunit + tSnap;
	while (timeunit < tStop)
	    {

	    // Evolve the LEBMs in the system for the minimum timestep in seconds
	    nBodySystem.evolveLEBMs(G,dtsec);

	    // Evolve the NBody particles for the minimum timestep in code units
	    nBodySystem.evolveSystem(dtunit);

	    timeunit = timeunit + dtunit;

	    // Recalculate the minimum timestep
	    dtunit = nBodySystem.calcCombinedTimestep();
	    dtsec = dtunit*unit2sec;

	    }

	printf("Time: %+.4E yr, Combined Timestep: %+.4E s, %+.4E units\n",timeunit/twopi, dtsec, dtunit);
	// Output data to files
	snapshotNumber++;
	timeyr = timeunit/twopi;

	// N Body data goes to a single file
	nBodySystem.outputNBodyData(outputfile, timeyr, input.orbitCentre);

	// LEBM data goes to separate files for each World in the System

	nBodySystem.outputLEBMData(snapshotNumber, timeyr, input.fullOutput);
	}

    //Close the file before returning

    fclose(outputfile);

    return 0;
    }

