# Written 28/3/14 by dh4gan
# Code reads in a sample of snapshots from EBM and plots them
# in a single image plot, where 
# y-axis = latitude
# x-axis = time

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


# Read log file to find time of all dumps

logfile = prefix+'.log'

time = np.genfromtxt(logfile,usecols=0)

# Select the time interval we are interested in
time = time[initialdump-1:finaldump]

print "Reading Files"

for i in range(nfiles):
    
    idump = initialdump + i    
    
    # Read in data
    inputfile = prefix+'.'+str(idump)
    outputfile = inputfile+'.png'
    
    data = np.genfromtxt(inputfile,skiprows=1)
    data = data[:,iy]
    
    #Obtain latitudes as well
    if(i==0):
        latitude = np.genfromtxt(inputfile,skiprows=1)
        latitude = latitude[:,1]
    
    # Add to complete dataset
    alldata.append(data)
    

print "File Read complete"

latitude = latitude*180.0/np.pi

# Convert this list of arrays into a single image

image = np.zeros((alldata[0].shape[0],nfiles))
                 
for i in range(len(alldata)):
    image[:,i] = alldata[i]
        

print "Plotting"

# Plot Image

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)', fontsize = 16)
ax.set_ylabel('Latitude (deg)', fontsize = 16)
ax.set_ylim(np.amin(latitude), np.amax(latitude))
ax.set_xlim(np.amin(time), np.amax(time))
plt.pcolor(time,latitude,image, cmap='spectral')
colourbar = plt.colorbar()
colourbar.set_label(namedict[ykey])


outputfile = prefix+'_'+ykey+'lat_vs_t.png'
plt.savefig(outputfile, format='png')
    