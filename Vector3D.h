/*
 * Vector3D.h
 *
 *  Created on: Mar 14, 2013
 *      Author: dhf
 *
 *      A simple 3D Vector Class (holds useful vector operations)
 */

using namespace std;

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include "math.h"
#include <vector>


class Vector3D{
    public:
	Vector3D();
	Vector3D(double x,double y,double z);
	virtual ~Vector3D();

	double dotProduct(Vector3D other);
	Vector3D scaleVector(double scale);
	double magVector();
	Vector3D unitVector();
	Vector3D addVector(Vector3D b);
	Vector3D addVector(Vector3D b, Vector3D c);
	Vector3D addVector(Vector3D b, Vector3D c, Vector3D d);
	Vector3D subtractVector(Vector3D b);
	Vector3D relativeVector(Vector3D b);
	Vector3D crossProduct(Vector3D other);
	Vector3D zeroVector();

	void rotateX(double angle);
	void rotateY(double angle);
	void rotateZ(double angle);

	void printVector();

	double elements [3];

};




#endif /* VECTOR3D_H_ */
