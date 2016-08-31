# Written 17/1/14 by dh4gan
# Code reads in a sample of snapshots from EBM and plots a time average
# latitudinal plot

import matplotlib.pyplot as plt
import numpy as np
from string import split
import io_oberon as io


# Find files to read in

prefix = raw_input("What is the name of the World? ")
t_interval = input ("What is the averaging interval (years)?")
dumpfreq = input("What is the snapshot interval (years)?")

initialdump = input("Which dump to start from? ")

# Read the header of the first file to find the dump frequency
inputfile = prefix+'.'+str(initialdump)
f = open(inputfile, 'r')
line = f.readline()

numbers = split(line)

#dumpfreq=float(numbers[1])       
#print dumpfreq, t_interval 

nfiles = np.int(t_interval/dumpfreq)+1
nfiles = 67
nlat = 144


print "Number of files is ", nfiles
# Define x axis as latitude (always)
xkey = 'lat'
ix = io.snapcoldict[xkey]

# Pick variable to time average
print "Which variable is to be time averaged?"

for i in range(len(io.snapvariablekeys)):
    if(i!=ix): print io.snapvariablekeys[i],":\t \t", io.snapnamedict[io.snapvariablekeys[i]]

keyword = raw_input("Enter appropriate keyword:   ")

ykey = keyword
iy = io.snapcoldict[keyword]

outputfile = ykey+'_timeaveraged_'+str(t_interval)+'yr.png'

alldata = []

for i in range(nfiles):
    
    idump = initialdump + i
    
    # Read in data
    inputfile = prefix+'.'+str(idump)
    
    data = np.genfromtxt(inputfile,skiprows=1)
    
    if(i==0):
        xdata = data[:,ix]
    
    data = data[:,iy]
    # Add to totals
    alldata.append(data)

# Calculate mean, standard deviation

mean= np.zeros((alldata[0].shape))
sd= np.zeros((alldata[0].shape))

nrows = alldata[0].shape[0]


for idump in range(nfiles):    
    mean[:] = mean[:] + alldata[idump][:]

mean = mean/nfiles

for idump in range(nfiles):    
        sd[:] = sd[:] + (mean[:]-alldata[idump][:])**2
        
sd = sd/(nfiles-1)
sd = np.sqrt(sd)

print "Mean and sd calculated"

# Plot analytical function

analytical = 302.3 - 45.3*np.sin(xdata)*np.sin(xdata)

# Make figure

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Latitude (degrees)', fontsize = 16)
ax.set_ylabel(io.snapnamedict[ykey], fontsize = 16)
ax.plot(xdata*180.0/np.pi,mean-1.0, label='Fiducial Model')
ax.plot(xdata*180.0/np.pi,analytical, linestyle='dashed', label='North & Coakley (1979) fit')
ax.legend(loc='lower center')
fig1.savefig(outputfile, format='png')
