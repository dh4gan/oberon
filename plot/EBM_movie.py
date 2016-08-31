# Written 17/1/14 by dh4gan
# Code reads in a sample of snapshots from EBM and plots them

import matplotlib.pyplot as plt
import numpy as np
from os import system

# File order
# 0 x
# 1 latitude
# 2 T
# 3 C
# 4 Q
# 5 IR
# 6 albedo
# 7 insolation
# 8 tau
# 9 ice fraction
# 10 habitability index


# Set up tuples and dictionaries

variablekeys = ("x","lat", "T", "C", "Q","IR","Albedo", "S", 
                "tau", "ice","hab")
variablenames = ("x", r"$\lambda$", "T (K)", "C (/)", r"Net Heating ($erg\,s^{-1}\, cm^{-2}$)",  
                 r"IR Cooling ($erg\,s^{-1}\, cm^{-2}$)",r"Albedo", r"Mean Insolation ($erg\,s^{-1}\, cm^{-2}$)", 
                 "Optical Depth", r"$f_{ice}$","Habitability Index")                  
variablecolumns = (0,1,2,3,4,5,6,7,8,9,10)

nvar = len(variablekeys)

namedict = {}
coldict = {}

for i in range(len(variablekeys)):
    namedict[variablekeys[i]] = variablenames[i]
    coldict[variablekeys[i]] = variablecolumns[i]


# Open log file and read contents

prefix = raw_input("What is the name of the World? ")


initialdump = input("Which dump to start from? ")
finaldump = input("Which is the final dump?")

nfiles = finaldump-initialdump +1

# Define x axis as latitude (always)
xkey = 'lat'
ix = coldict[xkey]

# Pick variable to time average
print "Which variable is to be plotted?"

for i in range(len(variablekeys)):
    if(i!=ix): print variablekeys[i],":\t \t", namedict[variablekeys[i]]

keyword = raw_input("Enter appropriate keyword:   ")

ykey = keyword
iy = coldict[keyword]

alldata = []

print "Plotting"
for i in range(nfiles):
    
    idump = initialdump + i
    
    print idump
    
    # Read in data
    inputfile = prefix+'.'+str(idump)
    outputfile = inputfile+'.png'
    
    data = np.genfromtxt(inputfile,skiprows=1)
    
    if(i==0):
        xdata = data[:,ix]
    
    data = data[:,iy]
    # Add to totals    


    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_xlabel(namedict[xkey], fontsize = 16)
    ax.set_ylabel(namedict[ykey], fontsize = 16)
    ax.plot(xdata,data)

    fig1.savefig(outputfile, format='ps')
    
# Make movie

# Command for converting images into gifs - machine dependent

#convertcommand = '/opt/ImageMagick/bin/convert '
convertcommand = '/usr/bin/convert '

print "Making movie"
system(convertcommand +prefix+'*.png movie.gif')
system('rm '+prefix+'*.png')
