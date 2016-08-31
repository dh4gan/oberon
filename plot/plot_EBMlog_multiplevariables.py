# Written 17/1/14 by dh4gan
# Code reads in log files from EBM and plots them

import matplotlib.pyplot as plt
import numpy as np
import filefinder.localfiles as ff
import io_oberon.io_EBM


# Find log files

inputfile = ff.find_local_input_files('*.log')
prefix = inputfile[:-4]

xkey, ykeys, ix, iylist, nyvar = io_oberon.io_EBM.select_multiple_variables(io_oberon.io_EBM.logvariablekeys, io_oberon.io_EBM.logcoldict, io_oberon.io_EBM.lognamedict)

# Read in data

print 'Reading ',inputfile

data = np.genfromtxt(inputfile)

# Find maximum and minimum time limits to plot

xmin = np.amin(data[:,ix])
xmax = np.amax(data[:,ix])

ymin = []
ymax = []
for iy in iylist:
    ymin.append(np.amin(data[:,iy]))
    ymax.append(np.amax(data[:,iy]))

ymin = np.array(ymin)
ymax = np.array(ymax)

print 'Data xmin, xmax: ', xmin, ' , ',xmax
for i in range(nyvar):
    print 'Data ymin, ymax: ', ymin[i], ' , ',ymax[i]

newlimits = raw_input("Set new x limits? (y/n) ")
if newlimits=='y':
    xmin = input("What is the minimum x value? ")
    xmax = input("What is the maximum x value? ")

    # Find x points in this range
    inxrange = np.where(np.logical_and(data[:,ix]>xmin,data[:,ix]<xmax))    
    
    # Find minimum and maximum y in this range
    
    for i in range(nyvar):
        ymin[i] = np.amin(data[inxrange,iylist[i]])
        ymax[i] = np.amax(data[inxrange,iylist[i]])
    
        print 'New '+ykeys[i]+' limits: ',ymin[i], ymax[i]

outputfile = ''

for i in range(nyvar):    
    outputfile = outputfile + ykeys[i]+'_'
    
outputfile = outputfile+'vs_'+xkey+'.png'

npoints = len(data[:,ix])

# If npoints large, prompt user for subsampling

nsubsample = 0
if(npoints > 1.0e4):
    nsubsample = input('Attempting to plot '+str(npoints)+' points: subsampling is recommended! \nEnter subsampling factor (enter 0 to skip):')
    nsubsample = int(nsubsample)

#globalymin = 0.0
#globalymax= 1.0
#if(ymax-ymin < 1.0):
#    ymax = ymax +0.5
#    ymin = ymin -0.5

# Make figure

print xmin, xmax

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel(io_oberon.io_EBM.lognamedict[xkey], fontsize = 16)
ax.set_ylabel(io_oberon.io_EBM.lognamedict[ykeys[0]], fontsize = 16)
ax.set_xlim(xmin,xmax)

if(nsubsample>0):
    lns1 = ax.plot(data[0::nsubsample,ix],data[0::nsubsample,iylist[0]],label=io_oberon.io_EBM.lognamedict[ykeys[0]], color='blue')
else:
    lns1 = ax.plot(data[:,ix],data[:,iylist[0]],label=io_oberon.io_EBM.lognamedict[ykeys[0]], color='blue')

ax2 = ax.twinx()

ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin[1],ymax[1])
ax2.set_ylabel(io_oberon.io_EBM.lognamedict[ykeys[1]], fontsize = 16)

if(nsubsample>0):
    lns2 = ax2.plot(data[0::nsubsample,ix],data[0::nsubsample,iylist[1]],label=io_oberon.io_EBM.lognamedict[ykeys[1]], color='green', linestyle='dashed')
else:
    lns2 = ax2.plot(data[:,ix],data[:,iylist[1]],label=io_oberon.io_EBM.lognamedict[ykeys[1]], color='green', linestyle='dashed')

lns = lns1+lns2
labs = [l.get_label() for l in lns]

ax2.legend(lns,labs,bbox_to_anchor=(0.,1.02,1.,1.02),loc=3,ncol=2,mode="expand",borderaxespad=0)
#ax2.legend(loc = 'upper right')

fig1.savefig(outputfile, format='png')

