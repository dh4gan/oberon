# Written 17/1/14 by dh4gan
# Code reads in log files from EBM and calculates Fourier Transforms

import matplotlib.pyplot as plt
import numpy as np
import io_oberon.io_EBM
import filefinder.localfiles as ff
from scipy.signal import blackman, periodogram


# Identify log file

inputfile = ff.find_local_input_files('*.log')
prefix = inputfile[:-4]

xkey = "t"
keyword = xkey
ix = io_oberon.io_EBM.logcoldict[keyword]

periodchoice = input("Enter a period you expect to see: ")

print "Which variable for the y-axis?"

for i in range(len(io_oberon.io_EBM.logvariablekeys)):
    if(i!=ix): print io_oberon.io_EBM.logvariablekeys[i],":\t \t", io_oberon.io_EBM.lognamedict[io_oberon.io_EBM.logvariablekeys[i]]

keyword = raw_input("Enter appropriate keyword:   ")

ykey = keyword
iy = io_oberon.io_EBM.logcoldict[keyword]

# Read in data

print 'Reading ',inputfile

data = np.genfromtxt(inputfile)

npoints = len(data[:,ix])

print 'Computing periodogram'

# Find maximum and minimum time limits to plot

xmin = np.amin(data[:,ix])
xmax = np.amax(data[:,ix])

# Compute the periodogram of y along the x axis    

dt = data[-1,ix]/npoints

# First, smooth and apply a window function

w = blackman(npoints)
smoothed = data[:,iy]

#fourier = fft(smoothed*w)
#print fourier
#freq = fftfreq(npoints,dt)
#fourier = fftshift(fourier)/npoints
#freq = fftshift(freq)
#powerspectrum = np.abs(fourier)**2

freqo, periodo = periodogram(smoothed*w, fs = 1.0/dt) 

ymin = np.amin(periodo)
ymax = np.amax(periodo)

# Calculate periods from these frequencies
periods = np.divide(1.0,freqo)

print 'Data xmin, xmax: ', xmin, ' , ',xmax
print 'Data ymin, ymax: ', ymin, ' , ',ymax

newlimits = raw_input("Set new x limits? (y/n) ")


if newlimits=='y':
    xmin = input("What is the minimum x value? ")
    xmax = input("What is the maximum x value? ")

    # Find x points in this range
    inxrange = np.where(np.logical_and(periods[:]>xmin,periods[:]<xmax))    
    
    # Find minimum and maximum y in this range
    
    ymin = np.amin(periodo[inxrange])
    ymax = np.amax(periodo[inxrange])
    
    print 'New y limits: ',ymin, ymax

xlines = np.zeros(6) 

for i in range(len(xlines)):
    xlines[i] = periodchoice/(i+1)

outputfile = 'fourier_'+ykey+'_vs_'+xkey+'.png'

# Make figure

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel("Period (years)", fontsize = 16)
ax.set_ylabel("Periodogram for "+io_oberon.io_EBM.lognamedict[ykey], fontsize = 16)

ax.plot(periods,periodo)

for i in range(len(xlines)):
    ax.axvline(x=xlines[i],color='r', linestyle = '--')

#plt.show()
fig1.savefig(outputfile, format='png')



