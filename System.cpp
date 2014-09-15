/*
 * System.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "System.h"
//#include "Constants.h"
#include <iostream>
#include <stdio.h>
#include <algorithm>

System::System()
    {
    name = "System";
    bodies = vector<Body*> (0);
    bodyCount = 0.0;
    Vector3D zeroVector;

    totalMass = 0.0;
    initialEnergy = 0.0;
    totalEnergy = 0.0;

    timeStep = 0.0;
    timeControl = 0.00002;

    G = 1.0;
    softeningLength = 1.0e-5;

    initialAngularMomentum = zeroVector;
    totalAngularMomentum = zeroVector;

    deltaAngularMomentum = 0.0;
    deltaEnergy = 0.0;

    positionCOM = zeroVector;
    velocityCOM = zeroVector;
    accelerationCOM = zeroVector;

    planetaryIlluminationOn = false;


    }

System::System(string &namestring, vector<Body*> &bodyarray)
    {

    name = namestring;
    bodies = bodyarray;
    bodyCount = bodies.size();
    Vector3D zeroVector;
    totalMass = 0.0;
    for (int i = 0; i < bodyCount; i++)
	{
	totalMass = totalMass + bodies[i]->getMass();
	}

    initialEnergy = 0.0;
    totalEnergy = 0.0;

    timeStep = 0.0;
    timeControl = 0.00002;

    G = 1.0;
    softeningLength = 1.0e-5;

    initialAngularMomentum = zeroVector;
    totalAngularMomentum = zeroVector;

    deltaAngularMomentum = 0.0;
    deltaEnergy = 0.0;

    positionCOM = zeroVector;
    velocityCOM = zeroVector;
    accelerationCOM = zeroVector;

    planetaryIlluminationOn = false;
    }

System::~System()
    {

    for (int i = 0; i << bodies.size(); i++)
	{
	delete bodies[i];
	bodies[i] = 0;
	}

    }

double safeAcos(double x) {
	if (x < -1.0)
		x = -1.0;
	else if (x > 1.0)
		x = 1.0;
	return acos(x);
}

/* Methods for Controlling Population of Bodies */

void System::addBody(Body* &newBody)
    {/* Author: dh4gan
     Adds a Body object to the System object */

    bodies.push_back(newBody);
    bodyCount = bodies.size();

    totalMass = totalMass + newBody->getMass();

    }

void System::removeBody(int bodyindex)
    {
    /* Author: dh4gan
     * removes Body (bodyindex) from the system
     */

    bodies.erase(bodies.begin() + bodyindex);

    // As a body has been deleted, must recalculate initial values of Energy, Angular Momentum
    calcInitialProperties();

    }


void System::setHostBodies(vector<int> orbitCentre)
    {
    /*
     * Written 19/8/14 by dh4gan
     * Sets up Host Bodies for Worlds where appropriate
     *
     */

    for (int i=0; i<bodyCount; i++)
	{
	if(bodies[i]->getType()=="World")
	    {
	    if(orbitCentre[i]>0)
		{
		bodies[i]->setHostBody(bodies[orbitCentre[i]-1]);
		}
	    }

	}



    }
/* Calculation Methods */

void System::calcCOMFrame(vector<int> participants)
    {/* Author: dh4gan
     * Calculates the position and velocity of the centre of mass
     * for a limited range of participants
     */

    int i;
    double m;
    Vector3D pos, vel, zerovector;

    positionCOM = zerovector;
    velocityCOM = zerovector;

    for (i = 0; i < bodyCount; i++)
	{
	if (participants[i] == 1)
	    {
	    pos = bodies[i]->getPosition();
	    vel = bodies[i]->getVelocity();
	    m = bodies[i]->getMass();

	    positionCOM = positionCOM.addVector(pos.scaleVector(m/totalMass));
	    velocityCOM = velocityCOM.addVector(vel.scaleVector(m/totalMass));

	    }
	}

    }



void System::calcCOMFrame()
    {/* Author: dh4gan
     Calculates the position and velocity of the centre of mass
     for all bodies */

    vector<int> participants(bodyCount,1);

    calcCOMFrame(participants);

    }
void System::transformToCOMFrame(vector<int> participants)
    {
    /* Author: dh4gan
     * Calculates Centre of Mass Frame and transforms the system to it
     * Calls the System method calcCOMFrame first */

    int i;
    Vector3D pos, vel;

    // Calculate Centre of Mass Vectors

    calcCOMFrame(participants);

    // Subtract these from Body vectors

    for (i = 0; i < bodyCount; i++)
	{
	if(participants[i]==1)
	    {
	pos = bodies[i]->getPosition();
	vel = bodies[i]->getVelocity();

	pos = pos.subtractVector(positionCOM);
	vel = vel.subtractVector(velocityCOM);

	bodies[i]->setPosition(pos);
	bodies[i]->setVelocity(vel);
	    }

	}
    //end of module
    }

void System::transformToCOMFrame()
    {
    /* Author: dh4gan
     * Calculates Centre of Mass Frame and transforms the system to it
     * Calls the System method calcCOMFrame first */

    vector<int> participants(bodyCount,1);

    transformToCOMFrame(participants);

    }


void System::transformToBodyFrame(int bodyIndex)
    {
    /* Author: dh4gan
     * Transforms system to the frame where body bodyIndex is at rest*/


    int i;
    Vector3D pos, vel;

    // Calculate Frame's vectors
    Vector3D framepos = bodies[bodyIndex]->getPosition();
    Vector3D framevel = bodies[bodyIndex]->getVelocity();


    // Subtract these from Body vectors

    for (i = 0; i < bodyCount; i++)
	{
	pos = bodies[i]->getPosition();
	vel = bodies[i]->getVelocity();

	pos = pos.subtractVector(framepos);
	vel = vel.subtractVector(framevel);

	bodies[i]->setPosition(pos);
	bodies[i]->setVelocity(vel);

	}
    //end of module
    }

void System::calcTotalEnergy()
    {
    /* Author: David the legend Harvey
     *
     * Purpose :
     * 		Calculate the total energy of the system
     *
     * Method :
     *      The total energy is the kinetic plus the gravitational potential
     *      so will need to loop through each body and work out the total of
     *      each
     *
     */

    int body;
    int other_body;
    int number_bodies;

    double r_distance;
    double body_ref_mass;
    double body_question_mass;
    double gravitational_potential = 0.0;
    double kinetic_energy = 0.0;
    double velocity;

    Vector3D r_vector;
    Vector3D body_ref_vel;
    Vector3D body_ref_pos;
    Vector3D body_question_pos;
    Vector3D r_hat;
    Vector3D force_total;

    //Should Test this as this could be faulty
    number_bodies = bodies.size();

    for (body = 0; body < number_bodies; body++)
	{

	//This is the body that is the reference for calc
	//the velocity
	body_ref_pos = bodies[body]->getPosition();
	body_ref_mass = bodies[body]->getMass();
	body_ref_vel = bodies[body]->getVelocity();

	for (other_body = 0; other_body < number_bodies; other_body++)
	    {
	    if (body != other_body)
		{
		//Need to loop through the other bodies in order to
		//find the gravitational potential

		body_question_pos = bodies[other_body]->getPosition();
		body_question_mass = bodies[other_body]->getMass();

		r_vector = body_question_pos.relativeVector(body_ref_pos);

		r_distance = r_vector.magVector();

		gravitational_potential += -G * body_ref_mass
			* body_question_mass / r_distance;

		}
	    }
	velocity = body_ref_vel.magVector();
	kinetic_energy += 0.5 * body_ref_mass * velocity * velocity;
	}

    totalEnergy = kinetic_energy + gravitational_potential;

    deltaEnergy = (totalEnergy - initialEnergy) / initialEnergy;

    }

void System::calcTotalAngularMomentum()
    {
    /*Author: dh4gan
     * Calculates the total Angular Momentum in the System, where the rotation axis is at the origin
     * Also calculates the deviation from the initial value.
     */

    Vector3D angularMomentum;
    Vector3D pos, vel, zerovector;

    double magTotalAngularMomentum, magInitialAngularMomentum;

    totalAngularMomentum = zerovector;

    for (int i = 0; i < bodyCount; i++)
	{
	pos = bodies[i]->getPosition().scaleVector(bodies[i]->getMass()); // multiply r by mass here
	vel = bodies[i]->getVelocity();

	totalAngularMomentum = totalAngularMomentum.addVector(pos.crossProduct(
		vel));

	}

    magTotalAngularMomentum = totalAngularMomentum.magVector();
    magInitialAngularMomentum = initialAngularMomentum.magVector();

    if (magInitialAngularMomentum != 0.0)
	{
	deltaAngularMomentum = (magInitialAngularMomentum
		- magTotalAngularMomentum) / magInitialAngularMomentum;
	}

    }

void System::calcNBodyTimestep(vector<Body*> &bodyarray, double dtmax)
    {/* Author: dh4gan
     Calls the calcTimestep method for every Body object in the System object,
     and finds the minimum value */

    int i;
    vector<double> dt(bodyCount,0.0);
    double dtarraymax;

#pragma omp parallel default(none) \
	shared(bodyarray,dt) \
	private(i)
	{
	for (i = 0; i < bodyCount; i++)
	    {

	    bodyarray[i]->calcTimestep(timeControl);

	    dt[i] = bodyarray[i]->getTimestep();
	    }

	}
    dtarraymax = *(min_element(dt.begin(), dt.end()));

    if(dtarraymax>dtmax)
	{
	dtarraymax = dtmax;
	}
    timeStep = dtarraymax;

    }


void System::setupOrbits(vector<int> orbitCentre)
    {

    /*
     * Written 22/1/14 by dh4gan
     * Given a list of orbit centres (i=0: CoM; i>0: body, i<0: (0,0,0))
     *
     */

    vector <int> participants(bodyCount,0);

    // Firstly, set up desired objects around centre of mass

    for (int b = 0; b < bodyCount; b++)
	{
	if (orbitCentre[b] == 0)
	    {
	    participants[b] = 1;

	    // Set up their orbits
	    bodies[b]->calcVectorFromOrbit(G,totalMass);
	    }
	if(orbitCentre[b]<0)
	    {

	    bodies[b]->calcVectorFromOrbit(G,totalMass-bodies[b]->getMass());
	    }


	}


    // Transform them to COM Frame
    transformToCOMFrame(participants);

    for (int b = 0; b < bodyCount; b++)
	{
	if (orbitCentre[b] > 0)
	    {

	    bodies[b]->calcVectorFromOrbit(G,
		    bodies[orbitCentre[b] - 1]->getMass());

	    Vector3D framepos =
		    bodies[orbitCentre[b] - 1]->getPosition().scaleVector(-1.0);
	    Vector3D framevel =
		    bodies[orbitCentre[b] - 1]->getVelocity().scaleVector(-1.0);
	    bodies[b]->changeFrame(framepos, framevel);


	    // Set up Hosts for Worlds (to calculate tidal heating)
	    if(bodies[b]->getType()=="World")
		{
		bodies[b]->setHostBody(bodies[orbitCentre[b]]);
		}
	    }

	}
    }


void System::calcForces(vector<Body*> &bodyarray)
    {
    /* Author: dh4gan
     This is a wrapper method, which abstracts the process of
     calculating accelerations, jerks, snaps and crackles of an array of bodies
     */

    int i;
    int length = bodyarray.size();
    Vector3D zeroVector;

#pragma omp parallel default(none) \
	shared(bodyarray,length,zeroVector) \
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < length; i++)
	    {
	    bodyarray[i]->setAcceleration(zeroVector);
	    bodyarray[i]->setJerk(zeroVector);
	    bodyarray[i]->calcAccelJerk(G, bodyarray, softeningLength);

	    }
	}
#pragma omp parallel default(none) \
	shared(bodyarray,length,zeroVector) \
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < length; i++)
	    {
	    bodyarray[i]->setSnap(zeroVector);
	    bodyarray[i]->setCrackle(zeroVector);
	    bodyarray[i]->calcSnapCrackle(G, bodyarray, softeningLength);
	    }
	}

    }

void System::calcInitialProperties()
    {
    /* Author: dh4gan
     * Calculates initial energy and angular momentum
     * and stores it for later evaluation
     *
     */

    // Update total number of bodies
    bodyCount = bodies.size();

    // Total Mass
    totalMass = 0.0;

    for (int i = 0; i < bodyCount; i++)
	{
	totalMass += bodies[i]->getMass();
	}

    // Put the system in the Centre of Mass Frame
    transformToCOMFrame();

    // Calculate initial forces
    calcForces(bodies);


    // Calculate Total Energy, and define initial value

    initialEnergy = totalEnergy;
    deltaEnergy = 0.0;


    //Calculate Total Angular Momentum, and define initial value
    calcTotalAngularMomentum();
    initialAngularMomentum = totalAngularMomentum;
    deltaAngularMomentum = 0.0;


    }

vector<double> System::checkForEclipses(int bodyIndex)
{
    /* Author: dh4gan
     This method looks from the position of the body given by bodyindex
     and checks whether any other bodies are obscured from view
     It does this by looking at impact parameters from the relative
     vector of a given body and bodyindex.  Any body which has an
     impact parameter smaller than its radius, and is in front of
     the given body, obscures it.

     The method returns a vector of doubles, describing what fraction of the star is eclipsed [0,1]
     */

	double pi = 3.14159265285;

	vector<double> eclipsefrac(bodyCount,0.0);
	Vector3D vector_i, vector_j;
	double mag_i, mag_j, idotj, b;
	double rad_i, rad_j, rad_i2, rad_j2;
	double angle1,angle2, area_i, area_j;

	for(int i=0; i<bodyCount; i++)
	{

		// Body cannot eclipse itself
		if(i==bodyIndex)
		{
			eclipsefrac[i]=0.0;
			continue;
		}

		//calculate relative vector between i and bodyindex

		vector_i = getBody(bodyIndex)->getPosition().relativeVector(getBody(i)->getPosition());
		mag_i = vector_i.magVector();

		// Get Radius of body i

		rad_i = getBody(i)->getRadius();

		// Now loop over other bodies to get impact parameters

		for (int j=0; j <bodyCount; j++)
		{
			// Skip for bodies i and bodyIndex

			if(j==i or j==bodyIndex)
			{
				continue;
			}

			// Calculate relative vector between body j and bodyIndex

			vector_j = getBody(bodyIndex)->getPosition().relativeVector(getBody(j)->getPosition());
			mag_j = vector_j.magVector();

			// Get Radius of body j

			rad_j = getBody(j)->getRadius();

			// impact parameter = mag(vector_j)sin alpha = mag(vector_j)*sqrt(1-(vector_i.vector_j)^2)

			if(mag_i*mag_j >0.0 and mag_i>mag_j)
			{
				idotj = vector_i.dotProduct(vector_j)/(mag_i*mag_j);
				b = mag_j*sqrt(1.0-idotj*idotj);


				// If impact parameter less than radius, and i further away than j, eclipse of i!
				// Calculate area covered during eclipse
				// (sum of two circular segment areas, one for each circle)

				if(b < (rad_i+rad_j) and idotj < 0.0)
				{

				    //b = b/rad_i;
				    //rad_i = 1.0;
				    //rad_j = rad_j/rad_i;


					rad_i2 = rad_i*rad_i;
					rad_j2 = rad_j*rad_j;
					angle1 = 2.0*safeAcos((rad_i2 + b*b -rad_j2)/(2.0*b*rad_i));
					angle2 = 2.0*safeAcos((rad_j2 + b*b -rad_i2)/(2.0*b*rad_j));

					area_i = 0.5*rad_i2*(angle1 - sin(angle1));
					area_j = 0.5*rad_j2*(angle2 - sin(angle2));

					eclipsefrac[i] = (area_i+area_j)/(pi*rad_i2);


					if(eclipsefrac[i]>1.0) {eclipsefrac[i] = 1.0;}
					if(eclipsefrac[i]<0.0) {eclipsefrac[i] = 0.0;}

				}
			}
			else
			{
				continue;
			}

		}// End loop over j to get impact parameters

	} //End loop over i to test eclipses


	return eclipsefrac;

}


void System::evolveSystem(double tbegin, double tend)
    {
    /* Author: dh4gan
     This method evolves the System of Body Objects using a 4th order Hermite integrator
     This is a predictor-corrector algorithm:
     the method stores predicted data in a vector of body objects
     */

    int i;
    double time;
    double dtmax = (tend-tbegin)/2.0;

    Vector3D pos, vel, acc, jerk, snap, crack; // Holders for the body state vectors
    Vector3D pos_p, vel_p, acc_p, jerk_p, snap_p, crack_p; // Holders for the predicted body state vectors
    Vector3D pos_c, vel_c, acc_c, jerk_c, snap_c, crack_c; // Holders for the predicted body state vectors

    Vector3D velterm, accterm, jerkterm;
    vector<Body*> predicted; // Vector to store the predicted body properties

    /* i. Calculate initial accelerations, jerks, snaps and crackles for all particles */

    calcInitialProperties();

    calcForces(bodies);

    /* ii. Calculate initial global timestep, total energy and total angular momentum */

    calcNBodyTimestep(bodies, dtmax);
    calcTotalEnergy();
    calcTotalAngularMomentum();

    /* iii Set predicted body array equal to the current body array */
    for (i=0;i<bodyCount; i++){
	predicted.push_back(bodies[i]->Clone());
    }

    /* Begin loop over time */
    time = tbegin;

    while (time < tend)
	{


#pragma omp parallel default(none) \
	shared(predicted)\
	private(i,pos,vel,acc,jerk,pos_p,vel_p)
	{
#pragma omp for schedule(runtime) ordered
	/* Calculate predicted positions and velocities */
	for (i = 0; i < bodyCount; i++)
	    {

	    // Pull the body object's data //
	    pos = bodies[i]->getPosition();
	    vel = bodies[i]->getVelocity();
	    acc = bodies[i]->getAcceleration();
	    jerk = bodies[i]->getJerk();


	    // 1. Calculate predicted position and velocity //
	    pos_p = pos.addVector(vel.scaleVector(timeStep), acc.scaleVector(
		    0.5 * timeStep * timeStep), jerk.scaleVector(timeStep
		    * timeStep * timeStep / 6.0));

	    vel_p = vel.addVector(acc.scaleVector(timeStep), jerk.scaleVector(
		    0.5 * timeStep * timeStep));

	    // Update the object holding the predicted data
	    predicted[i]->setPosition(pos_p);
	    predicted[i]->setVelocity(vel_p);

	    }
	}

	/* 2. Use predicted positions and velocities to calculate
	 * predicted accelerations, jerks, snaps and crackles */
	calcForces(predicted);



#pragma omp parallel default(none) \
	shared(predicted)\
	private(i,pos,vel,acc,jerk)\
	private(pos_p,vel_p,acc_p,jerk_p)\
	private(velterm,accterm,jerkterm,vel_c,pos_c)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < bodyCount; i++)
	    {

	    pos_p = predicted[i]->getPosition();
	    vel_p = predicted[i]->getVelocity();
	    acc_p = predicted[i]->getAcceleration();
	    jerk_p = predicted[i]->getJerk();

	    pos = bodies[i]->getPosition();
	    vel = bodies[i]->getVelocity();
	    acc = bodies[i]->getAcceleration();
	    jerk = bodies[i]->getJerk();

	    accterm = acc_p.addVector(acc).scaleVector(0.5 * timeStep);

	    jerkterm = jerk_p.relativeVector(jerk).scaleVector(timeStep
		    * timeStep / 12.0);

	    vel_c = vel.addVector(accterm, jerkterm);

	    accterm = acc_p.relativeVector(acc).scaleVector(timeStep * timeStep
		    / 12.0);
	    velterm = vel_c.addVector(vel).scaleVector(0.5 * timeStep);
	    pos_c = pos.addVector(velterm, accterm);

	    // update appropriate object with these new vectors

	    bodies[i]->setPosition(pos_c);
	    bodies[i]->setVelocity(vel_c);

	    }
	}

	/* 5. Calculate acceleration, jerk, snap and crackle for next step */

	calcForces(bodies);

	/* 5. Increase time by dt */
	time = time + timeStep;

	/* 6. Calculate new timestep, total energy, angular momentum and orbital data */
	calcNBodyTimestep(bodies, dtmax);
	calcTotalEnergy();
	calcTotalAngularMomentum();



	}
    // End of loop over time

    // Garbage collection on predicted pointers
    for(i=0; i< bodyCount; i++)
	{
	delete predicted[i];
	predicted[i]=0;
	}

    }

void System::evolveSystem(double dt)
    {

    double tbegin = 0.0;
    double tstop = dt;
    evolveSystem(tbegin,tstop);

    }

void System::calcPlanetaryEquilibriumTemperatures()
    {
    /*
     * Written 11/8/14 by dh4gan
     * Takes all Planets in the System, and calculates their equilibrium temperature
     * given the Stars in the System
     * Also calculates reflected starlight and thermal luminosity
     */

    double sep, temp, lum, rad, albedo;
    double pi = 3.141592653;
    double sigma_SB = 5.67e-8;

    double AU = 1.496e11;
    double lsol = 3.826e26;
    double rsol = 6.955e8;

    for (int j=0; j< bodyCount; j++)
	{
	if(bodies[j]->getType()=="Planet")
	    {

	    rad = bodies[j]->getRadius()*rsol/AU;
	    albedo = bodies[j]->getAlbedo();
	    temp = 0.0;
	    lum = 0.0;

	    for (int i=0; i<bodyCount; i++)
		{
		if(bodies[i]->getType() =="Star" and i!=j)
		    {

		    // Calculate separation between each star and planet
		    sep = bodies[j]->getPosition().subtractVector(bodies[i]->getPosition()).magVector();

		    // Add contribution to equilibrium temperature
		    temp = temp+ bodies[i]->getLuminosity()*lsol*(1.0-albedo)/(16.0*pi*sigma_SB*sep*sep*AU*AU);

		    //Calculate reflected starlight
		    lum = lum + bodies[i]->getLuminosity()*rad*rad*albedo/(sep*sep);

		    }
		}

	    temp = pow(temp,0.25);

	    // Set Planet's Equilibrium temperature
	    bodies[j]->setEquilibriumTemperature(temp);

	    // Set Planet's Reflective Luminosity
	    bodies[j]->setReflectiveLuminosity(lum);

	    // Calculate Planet's total luminosity
	    bodies[j]->calcLuminosity();

	    }

	}



    }

void System::evolveLEBMs(double &dt)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Checks for any World objects, and makes them update their LEBM models
     *
     */
    vector<double>eclipsefrac;

    if(planetaryIlluminationOn)
	{
	calcPlanetaryEquilibriumTemperatures();
	}

    for (int i=0; i< bodyCount; i++)
	{
	if(bodies[i]->getType()=="World")
	    {
	    eclipsefrac = checkForEclipses(i);
	    bodies[i]->updateLEBM(bodies,eclipsefrac,dt);

	    }

	}
    }



double System::calcCombinedTimestep()
    {
    /*
     * Written 10/1/14 by dh4gan
     * Calculates the minimum of the N Body and LEBM timesteps in the system
     *
     */

    double LEBMmin;
    double twopi = 2.0*3.141592;
    double tunit = 3.15e7/(twopi);

    double nbodymin = 1.0e30;
    // Calculate N Body minimum timestep in code units
    calcInitialProperties();
    calcForces(bodies);
    calcNBodyTimestep(bodies, nbodymin);

    nbodymin = getTimestep();

    // Calculate LEBM minimum timestep in seconds

    for (int b=0; b<bodyCount; b++)
	{

	if(bodies[b]->getType() == "World")
	    {
	    LEBMmin = bodies[b]->getLEBMTimestep();

	    LEBMmin =LEBMmin/tunit;

	    if(LEBMmin < nbodymin)
		{
		nbodymin = LEBMmin;
		}

	    }
	}
return nbodymin;

    }

void System::outputNBodyData(FILE *outputfile, double &time, vector<int> orbitCentre)
    {
    /*
     * Written 10/1/14 by dh4gan
     * Method writes N Body information to
     * already open file pointer
     * orbitCentre vector determines where the orbits are calculated from
     */

    Vector3D position, velocity;
    // Transform to the Centre of Mass Frame
    transformToCOMFrame();
    //transformToBodyFrame(0);

    for (int j = 0; j < bodyCount; j++)
	{

	if(orbitCentre[j]>0)
	    {
	    bodies[j]->calcOrbitFromVector(G, bodies[orbitCentre[j]-1]);
	    }

	else
	    {
	    bodies[j]->calcOrbitFromVector(G, totalMass);
	    }

	position = bodies[j]->getPosition();
	velocity = bodies[j]->getVelocity();

	//Write out the data
	//Output format  CSV
	// mass,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z

	fprintf(outputfile,
		"%+.4E,%+.4E, %s,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,"
			"%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E\n", time,
		totalEnergy, bodies[j]->getName().c_str(), bodies[j]->getMass(),
		bodies[j]->getRadius(), position.elements[0],
		position.elements[1], position.elements[2],
		velocity.elements[0], velocity.elements[1],
		velocity.elements[2], bodies[j]->getSemiMajorAxis(),
		bodies[j]->getEccentricity(), bodies[j]->getInclination(),
		bodies[j]->getLongitudeAscendingNode(),
		bodies[j]->getArgumentPeriapsis(), bodies[j]->getMeanAnomaly());

	}
    fflush(outputfile);

    }

void System::outputLEBMData(int &snapshotNumber, double &tSnap)
    {

    /*
     * Written 10/1/14 by dh4gan
     * Requires Worlds to write their LEBM data to files
     *
     */

    for (int b=0; b < bodyCount; b++)
	{
	if(bodies[b]->getType()=="World")
	    {

	    bodies[b]->outputLEBMData(snapshotNumber, tSnap);

	    }

	}

    }

