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
#include <chrono>
#include <time.h>
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
    printf("%s",screenBar.c_str());
    printf("\t\tOBERON - OBliquity and Energy balance Run on Nbody systems\n");
    printf("\t\tVersion: %s\n", VERSION);
    printf("\t\tDate: %s\n", __DATE__);
    printf("%s",screenBar.c_str());
    printf("  \n");
        
    // Record start time
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    

    // Read in parameters file

    if (argc == 2)
	{
	string fileString = string(argv[1]);
        
        printf("\tReading file %s\n", fileString.c_str());
        input.readFile(fileString);
	}
    else
	{
    
	input.readFile();
	
	}

        // Check and display input parameters
        input.checkParameters();
        input.displayParameters();

    // Record parameter data

    tMax = input.getDoubleVariable("MaximumTime");
    tSnap = input.getDoubleVariable("SnapshotTime");

        printf("Creating bodies \n");

    //First loop through each of the bodies in the system and set them up
    //adding them to the BodyArray
    for (i = 0; i < input.getIntVariable("Number_Bodies"); i++)
	{
        
        printf("Creating Body %s \n",input.getStringVariable("BodyName",i).c_str());
        if (input.getStringVariable("BodyType",i).compare("Star")==0)
        {
            BodyArray.push_back(new Star(input, i, G));
        }
        else if (input.getStringVariable("BodyType",i) == "Planet")
        {
        BodyArray.push_back(new Planet(input, i, G));
        }
        else if(input.getStringVariable("BodyType",i) == "World")
        {
            BodyArray.push_back(new World(input, i, G));
            
            if(input.getBoolVariable("Restart"))
            {
                cout << "Reading Temperature data for World " << BodyArray.back()->getName() << endl;
                snapshotNumber = BodyArray.back()->getRestartParameters();
        }
        
	    }

	printf("body %s set up \n", BodyArray.back()->getName().c_str());

	}

    // Set up System object using BodyArray
        
        string systemName = input.getStringVariable("SystemName");

    printf("Setting up system %s \n", systemName.c_str());

    nBodySystem = System(systemName, BodyArray);

    // If the System is created from orbital parameters, set up vectors here
        vector<int> orbitCentre;
        
        for (int i=0; i<input.getIntVariable("Number_Bodies"); i++)

        {
            orbitCentre.push_back(input.getIntVariable("OrbitCentre",i));
            
        }
    nBodySystem.setHostBodies(orbitCentre);

    if(input.getStringVariable("ParType").compare("Orbital")==0 and input.getBoolVariable("Restart")==false)
	{
	nBodySystem.setupOrbits(orbitCentre);
	}

    // Calculate the system's initial properties
    nBodySystem.calcInitialProperties();
    
    // Switch Planetary Illumination on/off
    nBodySystem.setIllumination(input.getBoolVariable("PlanetaryIllumination"));

    // Set up the outputs

    if (input.getBoolVariable("Restart") and snapshotNumber !=0)
	{
	outputfile = fopen(input.getStringVariable("NBodyOutput").c_str(), "a");
	}
    else
	{
	outputfile = fopen(input.getStringVariable("NBodyOutput").c_str(), "w");
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

    // This will be the default timestep measure - derive dtunit = dtsec*2pi/(3.15e7)
        
    dtsec = dtunit*unit2sec;
        
    timeunit = input.getDoubleVariable("SystemTime")*twopi; // system time in code units


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
	nBodySystem.outputNBodyData(outputfile, timeyr, orbitCentre);

	// LEBM data goes to separate files for each World in the System

	nBodySystem.outputLEBMData(snapshotNumber, timeyr, input.getBoolVariable("FullOutput"));
	}

    //Close the N Body file

    fclose(outputfile);

    
    // Simulation has now ended
    // Write elapsed runtime to the screen
    
    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
        
    printf("%s",screenBar.c_str());
    printf("Run %s complete \n", nBodySystem.getName().c_str());
    printf("Wall Clock Runtime: %f s \n", time_span.count());
    printf("%s",screenBar.c_str());
        
        
    return 0;
    }

