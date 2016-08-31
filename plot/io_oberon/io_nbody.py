# Written 15/1/14 Duncan Forgan
# Methods for reading nbody output

import numpy as np
import csv as c
import body as b

    
def read_nbody_datafile(filename, tmax):
    ''' Reads in nbody output to an array of body objects'''
    
    print "Reading ", filename        
    
    #Create the object file
    file_obj = c.reader( open(str(filename),"rb") )
    
    #First line should say the number of bodies separated by a comma
    number_bodies = np.int(file_obj.next()[1])
    
    #This is the format of the data file
    # tstop, etot, name, mass,radius, 
    # position_x,position_y,position_z,velocity_x,velocity_y,velocity_z
    # semimajor axis, eccentricity, inclination, longitude ascending node,
    # argument of periapsis, mean anomaly

    #Set up the bodies array
    bodyarray=[]
    
    #Initialise the bodies with initial parameters that it can be for n bodies
    for i in xrange(number_bodies):
        initial_params = file_obj.next()
        bodyarray.append(b.body(initial_params[2], initial_params[3], initial_params[4], \
                                   initial_params[5], initial_params[6], initial_params[7], \
                                   initial_params[8], initial_params[9], initial_params[10], \
                                   initial_params[11], initial_params[12], initial_params[13], \
                                   initial_params[14], initial_params[15], initial_params[16]))
    
    
    time=np.array(np.float(initial_params[0]))
    
    print "There are ",number_bodies," bodies in this system"
    if(tmax !=0.0): print "Plotting up to time ", tmax
    
    #Loop through the file and update the parameters of each of the bodies in question
    for row in file_obj:
        for i in xrange(number_bodies):
            if row[2] == bodyarray[i].bodytype:
                bodyarray[i].update_body(row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14],row[15],row[16])
        if row[2] == bodyarray[0].bodytype:
            time = np.append(time, np.float(row[0]))
        
        # Stop at the maximum time requested
        if time[-1]> tmax and tmax!=0.0:
            time = time[:-1] # delete last entry
            break
    
                
    print "There are ", np.size(time), " lines in the file"
    
    return time, bodyarray, number_bodies    


