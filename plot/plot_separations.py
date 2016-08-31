'''
Created on 7/3/14

@author: dh4gan

Show the positions of the bodies in the system

'''

from sys import argv
from matplotlib import pyplot as plt
import io_oberon.io_nbody
import numpy as np
    
# Data file can be read from the command line or from argument        
    
if len(argv)==1:
    input_file = raw_input("Enter the datafile: ")
    i1 = input("Enter a Number for Body 1: ")
    i2 = input("Enter a Number for Body 2: ")
elif len(argv)==2:
    input_file = argv[1]
    i1 = input("Enter a Number for Body 1: ")
    i2 = input("Enter a Number for Body 2: ")
elif len(argv)==3:    
    input_file = argv[1]
    i1 = int(argv[2])
    i2 = input("Enter a Number for Body 2: ")
elif len(argv)==4:
    input_file = argv[1]
    i1 = int(argv[2])
    i2 = int(argv[3])
    

    
tmax = 0.0
time, bodyarray, number_bodies = io_oberon.io_nbody.read_nbody_datafile(input_file, tmax)    
    
# Calculate separation of bodies 1 and 2

print bodyarray[i1].x

sepx = bodyarray[i2-1].x - bodyarray[i1-1].x 
sepy = bodyarray[i2-1].y - bodyarray[i1-1].y
sepz = bodyarray[i2-1].z - bodyarray[i1-1].z

sep = np.sqrt(sepx*sepx + sepy*sepy + sepz*sepz)

print sep

sepvx = bodyarray[i2-1].vx - bodyarray[i1-1].vx 
sepvy = bodyarray[i2-1].vy - bodyarray[i1-1].vy
sepvz = bodyarray[i2-1].vz - bodyarray[i1-1].vz

sepv = np.sqrt(sepvx*sepvx + sepvy*sepvy + sepvz*sepvz)


fig = plt.figure()
plt.subplot(311)
plt.xlabel("X Separation [ AU ]")
plt.ylabel("Y Separation [ AU ]")
plt.plot(sepx, sepy, '.', color='red')
        
plt.subplot(312)
plt.plot(time,sep,color='blue', label='$r$') 
plt.xlim(0,1000.0)   
plt.xlabel("Time [ years ]")
plt.ylabel("Separation [ AU ]")
plt.legend(loc='lower right')
    
plt.subplot(313)  
plt.xlim(0,1000.0)  
plt.plot(time, sepv,color='green', label='$v$')
plt.xlabel("Time [ years ]")
plt.ylabel("Relative Velocity [ AU / year ]")
plt.legend(loc='lower right')
    
    
plt.show()
        