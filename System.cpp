/*
 * System.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "System.h"
#include <iostream>
#include <stdio.h>

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
    positionCOM = zeroVector;
    velocityCOM = zeroVector;
    accelerationCOM = zeroVector;

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
    positionCOM = zeroVector;
    velocityCOM = zeroVector;
    accelerationCOM = zeroVector;

    }

System::~System()
    {

    for (int i = 0; i << bodies.size(); i++)
	{
	delete bodies[i];
	bodies[i] = 0;
	}

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

/* Calculation Methods */

void System::calcCOMFrame()
    {/* Author: dh4gan
     Calculates the position and velocity of the centre of mass */

    int i, j;
    double m;
    Vector3D pos, vel, zerovector;

    positionCOM = zerovector;
    velocityCOM = zerovector;

    for (i = 0; i < bodyCount; i++)
	{

	pos = bodies[i]->getPosition();
	vel = bodies[i]->getVelocity();
	m = bodies[i]->getMass();

	for (j = 0; j < 3; j++)
	    {
	    positionCOM.elements[j] = positionCOM.elements[j] + m
		    * pos.elements[j] / totalMass;
	    velocityCOM.elements[j] = velocityCOM.elements[j] + m
		    * vel.elements[j] / totalMass;
	    }
	}

    }
void System::transformToCOMFrame()
    {
    /* Author: dh4gan
     * Calculates Centre of Mass Frame and transforms the system to it
     * Calls the System method calcCOMFrame first */

    int i, j;
    Vector3D pos, vel;

    // Calculate Centre of Mass Vectors

    calcCOMFrame();

    // Subtract these from Body vectors

    for (i = 0; i < bodyCount; i++)
	{
	pos = bodies[i]->getPosition();
	vel = bodies[i]->getVelocity();

	for (j = 0; j < 3; j++)
	    {
	    pos.elements[j] = pos.elements[j] - positionCOM.elements[j];
	    vel.elements[j] = vel.elements[j] - velocityCOM.elements[j];

	    }
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

void System::calcGlobalTimestep(vector<Body*> &bodyarray, double dtmax)
    {/* Author: dh4gan
     Calls the calcTimestep method for every Body object in the System object,
     and finds the minimum value */

    int i;
    double dt;

    for (i = 0; i < bodyCount; i++)
	{
	bodyarray[i]->calcTimestep(timeControl);

	dt = bodyarray[i]->getTimestep();

	if (dt < dtmax and dt > 0.0)
	    {
	    dtmax = dt;
	    }
	}

    setTimestep(dtmax);

    }

void System::calcForces(vector<Body*> &bodyarray)
    {
    /* Author: dh4gan
     This is a wrapper method, which abstracts the process of
     calculating accelerations, jerks, snaps and crackles of an array of bodies
     */

    int length = bodyarray.size();
    Vector3D zeroVector;

    for (int i = 0; i < length; i++)
	{
	bodyarray[i]->setAcceleration(zeroVector);
	bodyarray[i]->setJerk(zeroVector);
	bodyarray[i]->calcAccelJerk(G, bodyarray, softeningLength);

	}

    for (int i = 0; i < length; i++)
	{
	bodyarray[i]->setSnap(zeroVector);
	bodyarray[i]->setCrackle(zeroVector);
	bodyarray[i]->calcSnapCrackle(G, bodyarray, softeningLength);
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
    calcTotalEnergy();
    initialEnergy = totalEnergy;
    deltaEnergy = 0.0;

    //Calculate Total Angular Momentum, and define initial value
    calcTotalAngularMomentum();
    initialAngularMomentum = totalAngularMomentum;
    deltaAngularMomentum = 0.0;

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

    calcGlobalTimestep(bodies, dtmax);
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

	/* 2. Use predicted positions and velocities to calculate
	 * predicted accelerations, jerks, snaps and crackles */
	calcForces(predicted);

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

	/* 5. Calculate acceleration, jerk, snap and crackle for next step */

	calcForces(bodies);

	/* 5. Increase time by dt */
	time = time + timeStep;

	/* 6. Calculate new timestep, total energy, angular momentum and orbital data */
	calcGlobalTimestep(bodies, dtmax);
	calcTotalEnergy();
	calcTotalAngularMomentum();

	for (i = 0; i < bodyCount; i++)
	    {
	    bodies[i]->calcOrbitFromVector(G, totalMass);
	    }

	}
    // End of loop over time

    // Garbage collection on predicted pointers
    for(i=0; i< bodyCount; i++)
	{
	delete predicted[i];
	predicted[i]=0;
	}

    }

