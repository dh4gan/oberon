# Written 17/1/14 by dh4gan
# Code reads in log files from multiple EBMs in the same N Body simulation and plots them

import matplotlib.pyplot as plt
import numpy as np
import io_oberon.io_EBM

import filefinder.localfiles as ff

# Find log files in directory

inputfiles = ff.find_sorted_local_input_fileset('*.log')

# Determine the plot of interest

xkey, ykey, ix,iy = io_oberon.io_EBM.select_variables(io_oberon.io_EBM.logvariablekeys, io_oberon.io_EBM.logcoldict, io_oberon.io_EBM.lognamedict)


# Prepare Figure

outputfile = ykey+'_vs_'+xkey+'.png'

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel(io_oberon.io_EBM.lognamedict[xkey], fontsize = 16)
ax.set_ylabel(io_oberon.io_EBM.lognamedict[ykey], fontsize = 16)

globalymin = 1.0e30
globalymax = -1.0e30

xlimcheck = False

for inputfile in inputfiles:
    print 'Reading ',inputfile

    # Determine the World Names by deleting the '.log' suffix
    worldname = inputfile.replace(".log", "")
    data = np.genfromtxt(inputfile)

    # Find maximum and minimum time limits to plot

    xmin = np.amin(data[:,ix])
    xmax = np.amax(data[:,ix])

    ymin = np.amin(data[:,iy])
    ymax = np.amax(data[:,iy])

    print 'Data xmin, xmax: ', xmin, ' , ',xmax
    print 'Data ymin, ymax: ', ymin, ' , ',ymax
    
                
    if(xlimcheck==False):
        newlimits = 'n'
        
        newlimits = raw_input("Set new x limits? (y/n) ")
        
    if newlimits=='y':
        if(xlimcheck==False):
            globalxmin = input("What is the minimum x value? ")
            globalxmax = input("What is the maximum x value? ")
        
        xlimcheck=True

        # Find x points in this range
        inxrange = np.where(np.logical_and(data[:,ix]>globalxmin,data[:,ix]<globalxmax))    
    
        # Find minimum and maximum y in this range
    
        ymin = np.amin(data[inxrange,iy])
        ymax = np.amax(data[inxrange,iy])
    
        print 'New y limits: ',ymin, ymax
    else:
        globalxmin= xmin
        globalxmax = xmax
    
    if(ymin < globalymin): globalymin = ymin    
    if(ymax > globalymax): globalymax = ymax
    
    ax.plot(data[:,ix],data[:,iy], label=worldname)

ax.set_xlim(globalxmin,globalxmax)
ax.set_ylim(globalymin,globalymax)
ax.legend()
fig1.savefig(outputfile, format='png')
