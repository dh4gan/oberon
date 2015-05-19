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

Input Options:


Outputs:

`<WorldName>.<number>` - a snapshot of <WorldName>'s latitudinal climate properties
`<WorldName>.log` - a log file for <WorldName> tracking global climate properties and position/orbital properties

An N Body file is also produced (name specified by the user), which is more suited to plotting the entire simulation's evolution (Star and Planet objects included).

Python scripts for plotting these datafiles can be found in +++

The code was developed using the eclipse CDT, which auto-generates a Makefile to compile the code.  There is also a manual Makefile in the repository to
compile with g++.
