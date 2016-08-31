'''
Created on Apr 16, 2013

@author: dh4gan

Plot the orbits of the bodies in the system

'''
import matplotlib
matplotlib.use('TkAgg')

from matplotlib import pyplot as plt
from sys import argv
import io_oberon.io_nbody


def plot_orbits():
    
    # Data file read from the command line        
    
    if len(argv)==1:
        input_file = raw_input("Enter the datafile: ")
    elif len(argv)==2:
        input_file = argv[1]
        tmax = input("What is the maximum time? Enter 0 for no maximum ")
    elif len(argv)==3:
        input_file = argv[1]
        tmax = argv[2]                                    
    
    time, bodyarray, number_bodies = io_oberon.io_nbody.read_nbody_datafile(input_file, tmax)
         
    # Now begin plotting data - add in an option to skip the star
    
    ignore_body1 = True
    
    if(len(argv)>3):
        if(argv[3]=='all'):
            ignore_body1 = False    
    
    # Plot r as a function of time
    
    rfig = plt.figure()
    rax = rfig.add_subplot(111)
    rax.set_xlabel("Time (yr)", fontsize= 16)
    rax.set_ylabel("r (AU)", fontsize= 16)
    
    # Plot semimajor axis as a function of time
    
    afig = plt.figure()
    aax = afig.add_subplot(111)
    aax.set_xlabel("Time (yr)", fontsize= 16)
    aax.set_ylabel("a (AU)", fontsize= 16)
    
    # Eccentricity
    
    efig = plt.figure()
    eax = efig.add_subplot(111)
    eax.set_xlabel("Time (yr)", fontsize= 16)
    eax.set_ylabel("Eccentricity", fontsize= 16)
    
    # Inclination
    
    ifig = plt.figure()
    iax = ifig.add_subplot(111)
    iax.set_xlabel("Time (yr)", fontsize= 16)
    iax.set_ylabel("Inclination", fontsize= 16)
    
    # Longitude of the Ascending Node
    
    longfig = plt.figure()
    longax = longfig.add_subplot(111)
    longax.set_xlabel("Time (yr)", fontsize= 16)
    longax.set_ylabel("Longitude of the Ascending Node", fontsize= 16)
    
    # Argument of Periapsis
    
    argfig = plt.figure()
    argax = argfig.add_subplot(111)
    argax.set_xlabel("Time (yr)", fontsize= 16)
    argax.set_ylabel("Argument of Periapsis", fontsize= 16)

    # Periapsis Radius
    
    rperfig = plt.figure()
    rperax = rperfig.add_subplot(111)
    rperax.set_xlabel("Time (yr)", fontsize= 16)
    rperax.set_ylabel("Periastron Radius (AU)", fontsize= 16)
    rperax.set_yscale('log')
        
    
    for i in xrange(number_bodies):
        if(i==0 and ignore_body1): continue # Skip star if necessary                

        rper = bodyarray[i].a*(1.0-bodyarray[i].e)
        
        rax.plot(time, bodyarray[i].r, label=bodyarray[i].bodytype)        
        aax.plot(time, bodyarray[i].a, label=bodyarray[i].bodytype)
        eax.plot(time, bodyarray[i].e, label=bodyarray[i].bodytype)
        iax.plot(time, bodyarray[i].i, label=bodyarray[i].bodytype)
        longax.plot(time, bodyarray[i].longascend, label=bodyarray[i].bodytype)
        argax.plot(time, bodyarray[i].argper, label=bodyarray[i].bodytype)
        rperax.plot(time,rper,label=bodyarray[i].bodytype)

    rax.legend(loc='upper left')
    aax.legend(loc='upper left')
    eax.legend(loc='upper left')
    iax.legend(loc='upper left')
    longax.legend(loc='upper left')
    argax.legend(loc='upper left')
    rperax.legend(loc='upper left')
    
    rfig.savefig('r.png', format='png')
    afig.savefig('a.png', format='png')
    efig.savefig('e.png', format='png')
    ifig.savefig('i.png', format='png')
    longfig.savefig('longascend.png', format='png')
    argfig.savefig('argper.png', format='png')
    rperfig.savefig('rperi.png', format='png')
    #plt.show()


if __name__ == '__main__':
    plot_orbits()
