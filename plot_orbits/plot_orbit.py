'''
Created on Apr 16, 2013

@author: David Harvey

Plot the orbits of the bodies in the system

'''
import matplotlib
matplotlib.use('TkAgg')

import body
import numpy as np
import csv as c
from matplotlib import pyplot as plt
import animate_orbit

def plot_orbits():
    
    #Where is the data file?
    input_dir = '/disk1/dhf/programs/cplusplus/nbody_EBM/Debug/'
    input_file = str(input_dir)+'nbody_output.txt'
    
    
    print "Reading ", input_file
    #Create the object file
    file_obj = c.reader( open(str(input_file),"rb") )
    
    #First line should say the number of bodies separated by a comma
    number_bodies = np.int(file_obj.next()[1])
    
    #This is the format of the data file
    # tstop, etot, name, mass,radius, 
    # position_x,position_y,position_z,velocity_x,velocity_y,velocity_z

    #Set up the bodies array
    bodyarray=[]
    
    #Initialise the bodies with initial params that it can be for n bodies
    for i in xrange(number_bodies):
        initial_params = file_obj.next()
        bodyarray.append(body.body(initial_params[2], initial_params[3], initial_params[4], \
                                   initial_params[5], initial_params[6], initial_params[7], \
                                   initial_params[7], initial_params[8], initial_params[9]))
    
    
    time=np.array(np.float(initial_params[0]))

    
    
    print "There are ",number_bodies," bodies in this system"
    #Loop through the file and update the paramters of each of the bodies in question
    for row in file_obj:
        for i in xrange(number_bodies):
            if row[2] == bodyarray[i].bodytype:
                bodyarray[i].update_body(row[5],row[6],row[7],row[8],row[9],row[10])
        if row[2] == bodyarray[0].bodytype:
            time = np.append(time, np.float(row[0]))

                
    print np.size( bodyarray[1].x)
    print np.size(time)
    
    for i in xrange(number_bodies):
        fig=plt.figure(i)
        plt.suptitle(str(bodyarray[i].bodytype))
        plt.subplot(211)
        plt.xlabel("X Position [ AU ]")
        plt.ylabel("Y Position [ AU ]")
        plt.plot(bodyarray[i].x, bodyarray[i].y, '.')
        
        plt.subplot(212)
        plt.plot(time, bodyarray[i].vy, '.')
        plt.xlabel("Time [ years ]")
        plt.ylabel("Velocity [ AU / year ]")

    plt.show()
        
    # First set up the figure, the axis, and the plot element we want to animate
    fig1 = plt.figure(number_bodies)
    
    plt.xlim(-1e-7, 1e-7 )
    plt.ylim(-1e-7, 1e-7)
    plt.xlabel('x')
    plt.title(str(bodyarray[0].bodytype))
    
    anim1 = animate_orbit.animate_orbit(bodyarray[0].x, bodyarray[0].y, fig1)
    plt.show()

    fig2 = plt.figure(number_bodies+1)

    plt.xlim(-1, 1 )
    plt.ylim(-1, 1)
    plt.xlabel('x')
    plt.title(str(bodyarray[1].bodytype))
    anim2 = animate_orbit.animate_orbit(bodyarray[1].x, bodyarray[1].y, fig2)
    
    plt.show()
    
        
    #Now loop through each line and 

if __name__ == '__main__':
    plot_orbits()
