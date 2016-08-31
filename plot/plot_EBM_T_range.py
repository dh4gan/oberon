# Written 17/1/14 by dh4gan
# Code reads in log files from EBM and plots them

import matplotlib.pyplot as plt
import numpy as np

import filefinder.localfiles as ff
from io_oberon.io_EBM import logcoldict, lognamedict

# Identify log file

inputfile = ff.find_local_input_files('*.log')
prefix = inputfile[:-4]

# Read in data

print 'Reading ',inputfile

data = np.genfromtxt(inputfile)

npoints = len(data[:,0])

# If npoints large, prompt user for subsampling

nsubsample = 0
if(npoints > 1.0e5):
    nsubsample = input('Attempting to plot '+str(npoints)+' points: subsampling is recommended! \nEnter subsampling factor (enter 0 to skip):')
    nsubsample = int(nsubsample)

print 'Plotting minimum, maximum and mean temperature for ',prefix

# Make three temperature figure

xkey = 't'
ix = logcoldict[xkey]
mincol = logcoldict['minT']
maxcol = logcoldict['maxT']
meancol = logcoldict['meanT']

xmin = np.amin(data[:,ix])
xmax = np.amax(data[:,ix])

inxrange = np.where(np.logical_and(data[:,ix]>xmin,data[:,ix]<xmax))

ymin = np.amin(data[inxrange,mincol])
ymax = np.amax(data[inxrange,maxcol])

print 'xmin, xmax: ', xmin, ' , ',xmax
print 'Global T min, max: ', ymin, ' , ',ymax

newlimits = raw_input("Set new x limits? (y/n) ")
if newlimits=='y':
    xmin = input("What is the minimum x value? ")
    xmax = input("What is the maximum x value? ")

    # Find x points in this range
    inxrange = np.where(np.logical_and(data[:,ix]>xmin,data[:,ix]<xmax))    
    
    # Find minimum and maximum y in this range
    
    ymin = np.amin(data[inxrange,mincol])
    ymax = np.amax(data[inxrange,maxcol])

    print 'New y limits: ',ymin, ymax



fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel(lognamedict[xkey], fontsize = 16)
ax.set_ylabel('T(K)', fontsize = 16)
ax.set_ylim(ymin,ymax)
ax.set_xlim(xmin,xmax)
if(nsubsample>0):
    ax.plot(data[1::nsubsample,ix],data[1::nsubsample,mincol], label='minimum T')
    ax.plot(data[1::nsubsample,ix],data[1::nsubsample,maxcol], label='maximum T')
    ax.plot(data[1::nsubsample,ix],data[1::nsubsample,meancol], label='mean T')

else:
    ax.plot(data[:,ix],data[:,mincol], label='minimum T')
    ax.plot(data[:,ix],data[:,maxcol], label='maximum T')
    ax.plot(data[:,ix],data[:,meancol], label='mean T')

ax.legend(loc='upper left')
fig1.savefig('compareT'+prefix+'.png', format='png')

