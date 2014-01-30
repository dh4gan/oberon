/*
 * Vector3D.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: dhf
 */


#include "Vector3D.h"
#include <iostream>

Vector3D::Vector3D()
    {
    elements[0] = elements[1] = elements[2] = 0.0;
    }

Vector3D::Vector3D(double x, double y, double z)
	    {
	elements[0] = x;
	elements[1] = y;
	elements[2] = z;
	    }


Vector3D::~Vector3D()
    {
    }

double Vector3D::dotProduct(Vector3D other)
// Written by dh4gan 11/3/13
// Simple method to take the dotProduct of two vectors
    {
    double dotProduct = 0.0;
    for (int i =0; i< 3; i++)
	{
	dotProduct += elements[i] *other.elements[i];
	}
    return dotProduct;
    }

Vector3D Vector3D::scaleVector(double scale)
// Written by dh4gan 11/3/13
// Simple method to scale a vector a by double b
    {
    Vector3D scaleVector;
    for (int i=0; i<3; i++)
	{
	scaleVector.elements[i] = elements[i] * scale;
	}

    return scaleVector;
    }

double Vector3D::magVector()
// Written by dh4gan 11/3/13
// Simple method to take the magnitude of a vector
    {

    double mag=0.0;

    for (int i=0; i<3; i++ )
	{
	mag += elements[i] * elements[i];
	}
    mag = sqrt(mag);
    return mag;
    }

Vector3D Vector3D::unitVector()
// Written by dh4gan 12/6/13
// Returns the unit vector
    {
    double mag = magVector();
    Vector3D unit;
    for (int i=0; i<3; i++)
	{
	unit.elements[i] = elements[i]/mag;
	}
    return unit;
    }

Vector3D Vector3D::addVector(Vector3D b)
// Written by dh4gan 11/3/13
// Simple method to add two vectors a and b
    {
    Vector3D add;

    for (int i=0;i < 3;i++)
	{
	add.elements[i] = elements[i] + b.elements[i];
	}

    return add;
    }

Vector3D Vector3D::addVector(Vector3D b, Vector3D c)
// Written by dh4gan 11/3/13
// Simple method to add two vectors a and b
    {
    Vector3D add;

    for (int i=0;i < 3;i++)
	{
	add.elements[i] = elements[i] + b.elements[i] + c.elements[i];
	}

    return add;
    }

Vector3D Vector3D::addVector(Vector3D b, Vector3D c, Vector3D d)
// Written by dh4gan 11/3/13
// Simple method to add two vectors a and b
    {
    Vector3D add;

    for (int i=0;i < 3;i++)
	{
	add.elements[i] = elements[i] + b.elements[i] + c.elements[i] + d.elements[i];
	}

    return add;
    }


Vector3D Vector3D::relativeVector(Vector3D b)
// Written by dh4gan 11/3/13
// Simple method to find the relative vector b-a
    {

    Vector3D relative;

    for (int i=0; i<3; i++)
	{
	relative.elements[i] = b.elements[i] - elements[i];
	}

    return relative;
    }

Vector3D Vector3D::subtractVector(Vector3D b)
// Written by dh4gan 11/3/13
// Simple method to find the relative vector b-a
    {

    Vector3D relative;

    for (int i=0; i<3; i++)
	{
	relative.elements[i] = elements[i] - b.elements[i];
	}

    return relative;
    }

Vector3D Vector3D::crossProduct(Vector3D other)
    {
    // Written by dh4gan 16/3/13
    // Calculates the cross product between a and b

    Vector3D cross;

    cross.elements[0] = elements[1] * other.elements[2] - elements[2] * other.elements[1];
    cross.elements[1] = elements[2] * other.elements[0] - elements[0] * other.elements[2];
    cross.elements[2] = elements[0] * other.elements[1] - elements[1] * other.elements[0];

    return cross;

    }


Vector3D Vector3D::zeroVector()
    {
    Vector3D zero;
    return zero;

    }

void Vector3D::rotateX(double angle)
    {
    // Written by dh4gan 1/8/13
    // Rotates a Vector around the x axis by angle

    Vector3D oldvec(elements[0],elements[1],elements[2]); // Use this variable to store previous position

    elements[0] = oldvec.elements[0];
    elements[1] = oldvec.elements[1]*cos(angle) - oldvec.elements[2]*sin(angle);
    elements[2] = oldvec.elements[1]*sin(angle) + oldvec.elements[2]*cos(angle);

    }

void Vector3D::rotateY(double angle)
    {
    // Written by dh4gan 1/8/13
    // Rotates a Vector around the y axis by angle

    Vector3D oldvec(elements[0],elements[1],elements[2]); // Use this variable to store previous position

    elements[0] = oldvec.elements[0]*cos(angle) + oldvec.elements[2]*sin(angle);
    elements[1] = oldvec.elements[1];
    elements[2] = -oldvec.elements[0]*sin(angle) + oldvec.elements[2]*cos(angle);

    }

void Vector3D::rotateZ(double angle)
    {
    // Written by dh4gan 1/8/13
    // Rotates a Vector around the z axis by angle

    Vector3D oldvec(elements[0],elements[1],elements[2]); // Use this variable to store previous position

    elements[0] = oldvec.elements[0]*cos(angle) - oldvec.elements[1]*sin(angle);
    elements[1] = oldvec.elements[0]*sin(angle) + oldvec.elements[1]*cos(angle);
    elements[2] = oldvec.elements[2];

    }

void Vector3D::printVector()
    {

    cout << elements[0] << "  " << elements[1] << "   " << elements[2] << endl;

    }
