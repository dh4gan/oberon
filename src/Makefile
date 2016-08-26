#####################################################
###                                               ###
### Makefile for oberon                           ###
###                                               ###
### Duncan H. Forgan 26/8/2014          	  ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables:
CC     = g++

# TODO - appropriate compiler flags?
CFLAGS = -O3

# Create object files:
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<
%.o: %.f90
	$(CC) $(CFLAGS) -c $<

# Source files (.F90)
SOURCESA = main.cpp Vector3D.cpp Body.cpp Star.cpp Planet.cpp World.cpp System.cpp parFile.cpp
OBJECTSA    = $(SOURCESA:.cpp=.o)

# Create executable files:
build: oberon

oberon:  $(OBJECTSA)
	$(CC) $(CFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o *.d oberon

# End Makefile
