# Written 17/1/14 by dh4gan
# Code reads in log files from multiple EBMs in the same simulation, and plots (meanT - mean(meanT))/mean(mean(T))

import matplotlib.pyplot as plt
import numpy as np
import io_oberon.io_EBM
import filefinder.localfiles as ff

# Open log files and read contents

inputfiles = ff.find_sorted_local_input_fileset('*.log')

# Select plotting variables

xkey, ykey, ix,iy = io_oberon.io_EBM.select_variables(io_oberon.io_EBM.logvariablekeys, io_oberon.io_EBM.logcoldict, io_oberon.io_EBM.lognamedict)


# Prepare Figure

outputfile = ykey+'_anomaly_vs_'+xkey+'.png'

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel(io_oberon.io_EBM.lognamedict[xkey], fontsize = 16)
ax.set_ylabel(io_oberon.io_EBM.lognamedict[ykey]+' Anomaly', fontsize = 16)


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

    # Calculate global mean of y data    
    ymean = np.mean(data[:,iy])
    anomaly = (data[:,iy]-ymean)
    
    # Calculate minimum and maximum anomaly    
    ymin = np.amin(anomaly)
    ymax = np.amax(anomaly)    

    print 'Data xmin, xmax: ', xmin, ' , ',xmax
    print 'Data ymin, ymax: ', ymin, ' , ',ymax
    print 'Mean y: ', ymean
    
                
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
            xdata = data[:,ix]
            
        # Calculate new mean of y data   
        ymean = np.mean(data[inxrange,iy])
        anomaly = data[:,iy]-ymean
        
        # Find minimum and maximum y in this range
    
        ymin = np.amin(anomaly[inxrange])
        ymax = np.amax(anomaly[inxrange])        
    
        print 'New ymin, ymax: ', ymin, ymax
        
    else:
        globalxmin= xmin
        globalxmax = xmax
    ax.plot(data[:,ix], anomaly, label=worldname)
    
    if(ymin < globalymin): globalymin = ymin    
    if(ymax > globalymax): globalymax = ymax        
    
print "Global x limits: ", globalxmin, globalxmax
print "Global y limits: ", globalymin, globalymax


ax.set_xlim(globalxmin,globalxmax)
ax.set_ylim(globalymin,globalymax)
ax.legend()

print 'Saving figure ', outputfile

fig1.savefig(outputfile, format='png')
