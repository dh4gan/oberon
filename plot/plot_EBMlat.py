# Written 17/12/15 by dh4gan
# Code reads in .lat file from nbody_EBM and plots T vs latitude vs time


import matplotlib.pyplot as plt
import numpy as np
import filefinder.localfiles as ff
import itertools


def moving_average(a,n=3):
    ret = np.cumsum(a,dtype=float)
    ret[n:] = ret[n:]-ret[:-n]
    return ret[n-1:]/n

# Open log file and read contents

inputfile = ff.find_local_input_files('*.lat')
prefix = inputfile[:-4]

# Read the header to get the latitudes

with open(inputfile) as t_in:
    latitude = np.genfromtxt(itertools.islice(t_in, 1))
    
latitude = latitude*180.0/np.pi

# Now read the rest of the file

print "Reading"
data = np.genfromtxt(inputfile, skip_header=1)
print "File Read"

# EDIT: Allow script to undersample by a factor of 10
# Do a moving average on the undersampled data

nsample = 10

# print 'Subsampling by a factor of ',nsample
# 
# fulldata = data
# 
# nlat = len(latitude)
# print nlat
# 
# for icol in range(nlat):
#     data[:,icol] = moving_average(fulldata[:,icol],n=nsample)
# 
# data = data[0:-1:10,:]


# Extract time and temperature data

time = data[:,0]
temperature = data[:,1:].transpose()

print "Plotting"

# Plot

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)', fontsize = 16)
ax.set_ylabel('Latitude (deg)', fontsize = 16)
ax.set_ylim(np.amin(latitude), np.amax(latitude))
ax.set_xlim(np.amin(time), np.amax(time))
plt.pcolor(time,latitude,temperature, cmap='spectral')
colourbar = plt.colorbar()
colourbar.set_label("Temperature (K)")

plt.savefig(prefix+'Tlat_vs_t.png')