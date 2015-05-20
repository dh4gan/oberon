This C++ code combines 1D latitudinal energy balance modelling of Earthlike climates with a
4th order Hermite N Body integrator.  

This code utilises object oriented programming - the N Body component of the code
uses the Body base class and the System class (which contains a vector of pointers to the Body objects in the simulation, and
methods for carrying out the N Body integration, and driver methods for carrying out the climate calculations).  The climate model component
of the code is encapsulated in three derived classes of Body: Star, Planet and  World (derived from Body).  The Star class provides thermal radiation
, and the Planet class provides reflected and thermal radiation based on the local radiation field.  Objects of the World class are considered to be Earthlike
in atmosphere and general composition, and their climates are evolved over the course of the simulation.

The code reads in a single input file (typically appended .params).  This parameter file contains 
a set of global parameters for all bodies in the simulation, along with specific parameters for each
body included in the simulation.

Parameter files can either specify the initial positions of all bodies, or the initial Keplerian orbits of all bodies.

The code is executed with the command

`> ./nbody_EBM input.params`


An example parameter file should be included with the repository

Input Options:

Global Options
--------------

ParType: 'Positional/Orbital' The format of data for entry (either Cartesian co-ordinates for position and velocity or orbital elements)
NBodyOutput: (string) The filename for the N Body data file
SnapshotTime: (double) The time interval between data dumps and snapshots
NGridPoints: (integer) The number of grid points used by the latitudinal energy balance model
MaximumTime: (double) The maximum runtime of the simulation
SystemName: (string) A descriptive string for the simulation
Number_Bodies: (integer) The number of bodies in the system
Restart: 'T/F' is the simulation a restart? True or False

Body Options (must be specified for all bodies)
------------
BodyName: (string)
BodyType: 'Star/Planet/World'
Mass: (double)
Radius: (double)

(Positional Files:
Position: (double) (double) (double) Cartesian position vector (AU) 
Velocity: (double) (double) (double) Cartesian velocity vector (2pi*AU/yr)

(Orbital Files:
SemiMajorAxis: (double) Semimajor Axis
Eccentricity: (double) Eccentricity
Inclination: (double) Inclination
LongAscend: (double) Longitude of the Ascending Node
Periapsis: (double) Argument of Periapsis
MeanAnomaly: (double) Mean Anomaly
OrbitCentre: (integer) Where is the initial orbit focus? -1 = (0,0,0), 0=system centre of mass, 1,2,3... = Body 1,2,3 )

Star Options
------------
Luminosity: (double) Bolometric Luminosity (solar luminosity)

World Options
-------------
RotationPeriod: (double) World rotation period in days
Obliquity: (double) World Obliquity in degrees
WinterSolstice: (double) orbital longitude of the winter solstice (degrees)
OceanFraction: (double, [0.0,1.0]) Fraction of the world's surface that is ocean
InitialTemperature: (double) Initial surface temperature of the world (at all latitudes)
IceMeltingOn: 'T/F' Is latent heat of melting for ice accounted for in climate calculation? True or False 


Outputs:
--------

`<WorldName>.<number>` - a snapshot of <WorldName>'s latitudinal climate properties
`<WorldName>.log` - a log file for <WorldName> tracking global climate properties and position/orbital properties

An N Body file is also produced (name specified by the user), which is more suited to plotting the entire simulation's evolution (Star and Planet objects included).

Python scripts for plotting these datafiles can be found in +++

The code was developed using the eclipse CDT, which auto-generates a Makefile to compile the code.  There is also a manual Makefile in the repository to
compile with g++.
