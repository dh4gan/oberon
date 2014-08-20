/*
 * Planet.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "Planet.h"

Planet::Planet() :
	Body()
    {
    albedo = 0.0;
    temperature = 0.0;
    reflectiveLuminosity = 0.0;
    luminosity = 0.0;
    }
Planet::Planet(string &namestring, string &typestring, double &m, double &rad,
	Vector3D &pos, Vector3D &vel, double &alb) :
	Body(namestring, typestring, m, rad, pos, vel)
    {
    albedo = alb;
    temperature = 0.0;
    reflectiveLuminosity = 0.0;
    luminosity = 0.0;
    }
Planet::Planet(string &namestring, string &typestring, double &m, double &rad,
	double semimaj, double ecc, double inc, double longascend,
	double argper, double meananom, double G, double totalMass, double &alb) :
	Body(namestring, typestring, m, rad, semimaj, ecc, inc, longascend,
		argper, meananom, G, totalMass)
    {
    albedo = alb;
    temperature = 0.0;
    reflectiveLuminosity = 0.0;
    luminosity = 0.0;
    }
Planet::~Planet()
    {
    }


void Planet::calcLuminosity()
    {
/*
 * Written 11/8/14 by dh4gan
 * Calculates the luminosity of the Planet, given its temperature and reflected starlight
 *
 */
    double pi = 3.141592653;
    double sigma_SB =5.67e-8;
    double AU = 1.496e11;
    double lsol = 3.8626e26;

    luminosity = 4.0*pi*getRadius()*getRadius()*AU*AU*sigma_SB*temperature*temperature*temperature*temperature/lsol;

    luminosity = luminosity + reflectiveLuminosity;




    }
