'''
Created on 7/3/14

@author: dh4gan

Show the positions of the bodies in the system

'''

from sys import argv
from matplotlib import pyplot as plt
import io_oberon.io_nbody
    
# Data file can be read from the command line or from argument        
    
if len(argv)==1:
    input_file = raw_input("Enter the datafile: ")
else:
    input_file = argv[1]
    
tmax = 0.0
time, bodyarray, number_bodies = io_oberon.io_nbody.read_nbody_datafile(input_file, tmax)    
    
for i in xrange(number_bodies):
    fig = plt.figure(i)
    plt.suptitle(str(bodyarray[i].bodytype))
    plt.subplot(211)
    plt.xlabel("X Position [ AU ]")
    plt.ylabel("Y Position [ AU ]")
    plt.plot(bodyarray[i].x, bodyarray[i].y, '.', color='red')
        
    plt.subplot(212)
    plt.plot(time,bodyarray[i].vx, '.',color='blue', label='$v_x$')
    plt.plot(time, bodyarray[i].vy, '.',color='green', label='$v_y$')
    plt.xlabel("Time [ years ]")
    plt.ylabel("Velocity [ AU / year ]")
    plt.legend(loc='lower right')
    
    
plt.show()
        