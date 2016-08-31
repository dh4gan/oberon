# Written 15/1/14 Duncan Forgan
# Useful methods for writing nbody_cplusplus parameter files

import numpy as np
import csv as c
import body as b
#from IPython.core.release import keywords

# EBM Log File dictionaries

logvariablekeys = ("t","minT", "maxT", "meanT", "meanQ","meanA","MeanIR", "meanS", 
                "meanhab", "meantidal", "a","ecc","inc","argper","longascend", "meananom", "obliq", "precession", "co2", "diffusion")
logvariablenames = ("Time (yr)", "minimum T (K)", "maximum T (K)", "mean T (K)", r"Net Heating ($erg\,s^{-1}\, cm^{-2}$)", r"Albedo", 
                 r"Mean Cooling ($erg\,s^{-1}\, cm^{-2}$)", r"Mean Insolation ($erg\,s^{-1}\, cm^{-2}$)", "Habitability Index", "Tidal Heating",
                 "Semimajor Axis", "Eccentricity", "Inclination","Argument of Periapsis", "Longitude of the Ascending Node",  "Mean Anomaly", "Obliquity", "Precession", r"$P_{CO_2}$", r"Diffusion Constant")
logvariablecolumns = range(20)


lognamedict = dict(zip(logvariablekeys,logvariablenames))
logcoldict = dict(zip(logvariablekeys,logvariablecolumns))


# EBM snapshot dictionaries

snapvariablekeys = ("x","lat", "T", "C", "Q","IR","Albedo", "S", 
                "tau", "ice","hab")
snapvariablenames = ("x", r"$\lambda$", "T (K)", "C (/)", r"Net Heating ($erg\,s^{-1}\, cm^{-2}$)",  
                 r"IR Cooling ($erg\,s^{-1}\, cm^{-2}$)",r"Albedo", r"Mean Insolation ($erg\,s^{-1}\, cm^{-2}$)", 
                 "Optical Depth", r"$f_{ice}$","Habitability Index")                  
snapvariablecolumns = (0,1,2,3,4,5,6,7,8,9,10)

snapnamedict = dict(zip(snapvariablekeys,snapvariablenames))
snapcoldict = dict(zip(snapvariablekeys,snapvariablecolumns))


# Methods for reading N body files

def checkData(N,names, mass, radius, semimaj, eccentricity, inclination, longascend, periapsis, meananomaly ):
    if len(names)!=N: 
        print "Error! Not enough names"
        return -1
    
    if len(mass)!=N: 
        print "Error! Not enough mass"
        return -1
    
    if len(radius)!=N: 
        print "Error! Not enough radius"
        return -1
    
    if len(semimaj)!=N: 
        print "Error! Not enough semimajor axes"
        return -1
    
    if len(eccentricity)!=N: 
        print "Error! Not enough eccentricities"
        return -1
    
    if len(inclination)!=N: 
        print "Error! Not enough inclinations"
        return -1
    
    if len(longascend)!=N: 
        print "Error! Not enough longitudes of ascending node"
        return -1
    
    if len(periapsis)!=N: 
        print "Error! Not enough periapses"
        return -1
    
    if len(meananomaly)!=N: 
        print "Error! Not enough mean anomalies"
        return -1
    
    return 0

def writeOrbitalParameterFile(paramfile, outputfile, tmax, tsnap, systemname, N, names, mass, radius, semimaj, eccentricity, inclination, longascend, periapsis, meananomaly):
    '''Writes parameter files which defines the initial orbits of objects'''
    
    # Firstly, check that each array has the correct number of entries
    
    success = checkData(N, names, mass, radius, semimaj, eccentricity, inclination, longascend, periapsis, meananomaly)
    
    if(success !=0):
        print "Failure in file write to "+paramfile
        return
    
    f = open(paramfile, 'w')
    
    line = "ParType \t Orbital \n"
    f.write(line)
    
    line = "OutputFile \t "+outputfile+ " \n"
    f.write(line)
        
    line = "SystemName \t "+systemname + " \n"
    f.write(line)
    
    line = "Number_Bodies \t "+str(N) + " \n"
    f.write(line)
    
    line = "MaximumTime \t "+str(tmax)+ " \n"
    f.write(line)
    
    line = "SnapshotTime \t "+str(tsnap)+ " \n"
    f.write(line)
    
    for i in range(N):
    
        line = "BodyType \t "+names[i] + " \n"
        f.write(line)
        
        line = "Mass \t "+str(mass[i]) + " \n"
        f.write(line)
        
        line = "Radius \t "+str(radius[i]) + " \n"
        f.write(line)
        
        line = "SemiMajorAxis \t "+str(semimaj[i]) + " \n"
        f.write(line)
        
        line = "Eccentricity \t "+str(eccentricity[i]) + " \n"
        f.write(line)
        
        line = "Inclination \t "+str(inclination[i]) + " \n"
        f.write(line)
        
        line = "LongAscend \t "+str(longascend[i]) + " \n"
        f.write(line)
        
        line = "Periapsis \t "+str(periapsis[i]) + " \n"
        f.write(line)
        
        line = "MeanAnomaly \t "+str(meananomaly[i]) + " \n"
        f.write(line)
    f.close()
    
def read_nbody_datafile(filename, tmax):
    ''' Reads in nbodycplusplus output to an array of body objects'''
    
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


def select_variables(variablekeys, coldict, labeldict):
    print "Which variable for the x-axis?"

    for i in range(len(variablekeys)):
        print variablekeys[i],":\t \t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    xkey = keyword
    ix = coldict[keyword]

    print "Which variable for the y-axis?"

    for i in range(len(variablekeys)):
        if(i!=ix): print variablekeys[i],":\t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    ykey = keyword
    iy = coldict[keyword]
        
    return xkey, ykey, ix,iy


def select_multiple_variables(variablekeys, coldict, labeldict):
    print "Which variable for the x-axis?"

    for i in range(len(variablekeys)):
        print variablekeys[i],":\t \t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    xkey = keyword
    ix = coldict[keyword]

    print "Which variables for the y-axis?"

    for i in range(len(variablekeys)):
        if(i!=ix): print variablekeys[i],":\t \t", labeldict[variablekeys[i]]

    keywords = raw_input("Enter keywords separated by spaces:   ")

    # Otherwise, parse keyword string into individual choices
    ykeys = keywords.split()
        
    
    iylist = []
    for word in ykeys:
        iylist.append(coldict[word])
    
    nyvar = len(iylist)        
    return xkey, ykeys, ix, iylist,nyvar


