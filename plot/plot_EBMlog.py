# Written 17/1/14 by dh4gan
# Code reads in log files from EBM and plots them

import matplotlib.pyplot as plt
import numpy as np

import filefinder.localfiles as ff
import io_oberon.io_EBM

# Identify log file

inputfile = ff.find_local_input_files('*.log')
prefix = inputfile[:-4]

# Select plotting variables

xkey, ykey, ix,iy = io_oberon.io_EBM.select_variables(io_oberon.io_EBM.logvariablekeys, io_oberon.io_EBM.logcoldict, io_oberon.io_EBM.lognamedict)

# Read in data

print 'Reading ',inputfile

data = np.genfromtxt(inputfile)

# Find maximum and minimum time limits to plot

xmin = np.amin(data[:,ix])
xmax = np.amax(data[:,ix])

ymin = np.amin(data[:,iy])
ymax = np.amax(data[:,iy])


print 'Data xmin, xmax: ', xmin, ' , ',xmax
print 'Data ymin, ymax: ', ymin, ' , ',ymax

newlimits = raw_input("Set new x limits? (y/n) ")
if newlimits=='y':
    xmin = input("What is the minimum x value? ")
    xmax = input("What is the maximum x value? ")

    # Find x points in this range
    inxrange = np.where(np.logical_and(data[:,ix]>xmin,data[:,ix]<xmax))    
    
    # Find minimum and maximum y in this range
    
    ymin = np.amin(data[inxrange,iy])
    ymax = np.amax(data[inxrange,iy])
    
    print 'New y limits: ',ymin, ymax


outputfile = ykey+'_vs_'+xkey+'.png'


npoints = len(data[:,iy])


# If npoints large, prompt user for subsampling

nsubsample = 0
if(npoints > 1.0e5):
    nsubsample = input('Attempting to plot '+str(npoints)+' points: subsampling is recommended! \nEnter subsampling factor (enter 0 to skip):')
    nsubsample = int(nsubsample)

# Make figure

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel(io_oberon.io_EBM.lognamedict[xkey], fontsize = 16)
ax.set_ylabel(io_oberon.io_EBM.lognamedict[ykey], fontsize = 16)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
if(nsubsample>0): 
    ax.plot(data[1::nsubsample,ix],data[1::nsubsample,iy])
else:
    ax.plot(data[:,ix],data[:,iy])

fig1.savefig(outputfile, format='png')

