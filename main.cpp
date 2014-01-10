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

    double G = 1;
    double year = 3.15e7;
    double pi = 3.141592654;
    double twopi = 2.0*pi;
    int i, j, fileType;

    double tBegin, tStop;
    double timeunit;
    double dtsec,dtunit;
    double tSnap, tMax;
    double totalEnergy;



    string body_i_name;
    double body_i_mass;
    double body_i_radius;
    Vector3D body_i_position;
    Vector3D body_i_velocity;

    System nBodySystem;
    vector<Body*> BodyArray;
    vector<Body*> Bodies;

    parFile input;
    FILE * outputfile, *accelfile;

    printf("  \n");
    printf("*********************************************** \n");
    printf("    NBODY EBM Code \n ");
    printf("    Date Created : 9th January 2014 \n");
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

    // Record parameter data

    tMax = input.maximumTime;
    tSnap = input.outputFrequency;

    //First loop through each of the bodies in the system and set them up
    //adding them to the BodyArray
    for (i = 0; i < input.number_bodies; i++)
	{
	if (fileType == 0)
	    {

	    cout << "setting up body from vectors" << endl;

	    body_i_position = input.getBodyPosition(i);
	    body_i_velocity = input.getBodyVelocity(i);

	    // If the Body is a Star, add a Star Object

	    if (input.BodyTypes[i] == "Star")
		{
		BodyArray.push_back(
			new Star(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.luminosity[i]));
		}

	    // If the Body is a Planet, add a Planet Object
	    if (input.BodyTypes[i] == "Planet")
		{
		BodyArray.push_back(
			new Planet(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity));
		}

	    // If the Body is a World, add a World Object and set up LEBM

	    if (input.BodyTypes[i] == "World")
		{
		BodyArray.push_back(
			new World(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.nPoints, input.obliquity[i], input.winterSolstice[i],
				input.oceanFraction[i], input.initialTemperature[i]));
		}


	    }
	else if (fileType == 1)
	    {
	    cout << "setting up body with orbital parameters" << endl;

	    // If the Body is a Star, add a Star Object
	    if (input.BodyTypes[i] == "Star")
		{
		BodyArray.push_back(
			new Star(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass, input.luminosity[i]));

		}

	    // If the Body is a Planet, add a Planet Object
	    if (input.BodyTypes[i] == "Planet")
		{
		BodyArray.push_back(
			new Planet(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass));

		}

	    // If the Body is a World, add a World Object and set up LEBM
	    if (input.BodyTypes[i] == "World")
		{
		BodyArray.push_back(
			new World(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass,input.nPoints, input.obliquity[i], input.winterSolstice[i],
				input.oceanFraction[i], input.initialTemperature[i]));
		}

	    }

	cout << "body " << BodyArray.back()->getName() << " set up" << endl;

	}

    // Set up System object using BodyArray

    cout << "Setting up system " << input.SystemName << endl;
    nBodySystem = System(input.SystemName, BodyArray);

    // Set up the outputs

    outputfile = fopen("nbody_output.txt", "w");
    fprintf(outputfile, "Number of Bodies, %i \n", input.number_bodies);

    accelfile = fopen("nbody_acc.txt", "w");

    // Now loop over snap shots, outputting the system data each time

    tStop = 0.0;
    tMax = tMax * 2.0 * pi;

    // Timesteps will be calculated in NBody units, and converted back to seconds for LEBM

    // Calculate the minimum LEBM timestep for all worlds and NBody timestep
    dtunit = nBodySystem.calcCombinedTimestep();

    dtsec = dtunit*twopi*3.15e7; // This will be the default timestep measure - derive dtunit = dtsec*2pi/(3.15e7)


    while (timeunit < tMax)
	{
	tBegin = tStop;
	tStop = timeunit + tSnap;

	while (timeunit < tStop)
	    {
	    // Evolve the LEBMs in the system for the minimum timestep in seconds
	    nBodySystem.evolveLEBMs(dtsec);

	    // Evolve the NBody particles for the minimum timestep in code units
	    nBodySystem.evolveSystem(dtunit);

	    // Recalculate the minimum timestep
	    dtunit = nBodySystem.calcCombinedTimestep();
	    dtsec = dtunit*twopi*3.15e7;
	    }

	// Output files TODO

	// Transform to the Centre of Mass Frame
	nBodySystem.transformToCOMFrame();
	Bodies = nBodySystem.getBodies();

	totalEnergy = nBodySystem.getTotalEnergy();

	for (j = 0; j < input.number_bodies; j++)
	    {

	    body_i_name = Bodies[j]->getName();
	    body_i_mass = Bodies[j]->getMass();
	    body_i_radius = Bodies[j]->getRadius();

	    body_i_position = Bodies[j]->getPosition();
	    body_i_velocity = Bodies[j]->getVelocity();

	    //Write out the data
	    //Output format  CSV
	    // mass,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z

	    fprintf(outputfile,
		    "%+.4E,%+.4E, %s,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E\n",
		    tStop, totalEnergy, body_i_name.c_str(), body_i_mass,
		    body_i_radius, body_i_position.elements[0],
		    body_i_position.elements[1], body_i_position.elements[2],
		    body_i_velocity.elements[0], body_i_velocity.elements[1],
		    body_i_velocity.elements[2]);
	    fflush(outputfile);
	    }

	fprintf(accelfile, "%+.4E %+.4E %+.4E %+.4E %+.4E\n", tStop,
		nBodySystem.getTimestep(), totalEnergy,
		Bodies[1]->getAcceleration().magVector(),
		Bodies[1]->getJerk().magVector());
	fflush(accelfile);
	}

    //Close the file before returning

    fclose(outputfile);
    fflush(accelfile);

    return 0;
    }

