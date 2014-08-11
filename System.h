/*
 * System.h
 *
 *  Created on: Nov 8, 2012
 *      Author: dh4gan
 *
 *      This is the header file setting up the System Class
 */

using namespace std;

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "Body.h"


class System {
public:
	System();
	System(string &namestring, vector<Body*> &bodyarray);
	virtual ~System();
	/* Accessor methods */

	string getName() {return name;}
	int getBodyCount(){ return bodyCount;}
	double getTotalMass() { return totalMass;}
	double getTimestep() { return timeStep;}
	double getTotalEnergy() { return totalEnergy; }
	vector<Body*> getBodies(){return bodies;}
	Body* getBody(int index) {return bodies[index];}
	Vector3D getpositionCOM(){ return positionCOM;}
	Vector3D getvelocityCOM(){ return velocityCOM;}
	Vector3D getaccelerationCOM(){ return accelerationCOM;}


	/* Variable Setting Methods */

	void setName(string namestring){name=namestring;}
	void setBodyCount(int count){ bodyCount=count;}
	void setTotalMass(double mtot) {totalMass = mtot;}
	void setTimestep(double dt) {timeStep = dt;}

	void setBodies(vector<Body*> bod){bodies=bod;}

	void setpositionCOM(Vector3D  r){positionCOM=r;}
	void setvelocityCOM(Vector3D  v){ velocityCOM=v;}
	void setaccelerationCOM(Vector3D  a){ accelerationCOM=a;}

	// Standard cloning method
	virtual System* Clone() { return new System(*this); }

	// Calculation Methods

	void addBody(Body* &newBody);
	void removeBody(int bodyindex);

	void calcNBodyTimestep(vector<Body*> &bodyarray, double dtmax);
	double calcCombinedTimestep();

	void calcCOMFrame(vector<int> participants);
	void transformToCOMFrame(vector<int> participants);
	void calcCOMFrame();
	void transformToCOMFrame();
	void transformToBodyFrame(int bodyIndex);

	void setupOrbits(vector<int> bodyCentre);

	void calcTotalEnergy();
	void calcTotalAngularMomentum();
	void calcInitialProperties();
	void calcForces(vector<Body*> &bodyarray);
	vector<double> checkForEclipses(int bodyindex);
	void evolveSystem(double tbegin, double tend);
	void evolveSystem(double dt);

	void calcPlanetaryEquilibriumTemperatures();
	void evolveLEBMs(double &dt);

	// Output Methods

	void outputNBodyData(FILE* outputfile, double &time, vector<int>orbitCentre);
	void outputLEBMData(int &snapshotNumber, double &tSnap);

protected:

	// Basic variables

	string name;
	vector<Body*> bodies;

	int bodyCount;

	double totalMass;
	double initialEnergy;
	double totalEnergy;

	double timeStep;
	double timeControl;

	double G;
	double softeningLength;

	double deltaAngularMomentum;  // Change in Angular Momentum since simulation begins
	double deltaEnergy; // Change in Energy since simulation begins

	Vector3D initialAngularMomentum;
	Vector3D totalAngularMomentum;

	Vector3D positionCOM;  // Position of the Centre of Mass
	Vector3D velocityCOM;  // Velocity of the Centre of Mass
	Vector3D accelerationCOM; // Acceleration of the Centre of Mass


};


#endif /* SYSTEM_H_ */
