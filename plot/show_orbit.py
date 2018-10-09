'''
Created on Apr 16, 2013

@author: David Harvey

Plot the orbits of the bodies in the system

'''
import matplotlib
matplotlib.use('TkAgg')


import numpy as np
from sys import argv
from matplotlib import pyplot as plt
import io_oberon.io_nbody
from time import sleep

def show_orbit():
    
    # Data file read from the command line        
    
    if len(argv)==1:
        input_file = raw_input("Enter the datafile: ")
    else:
        input_file = argv[1]
    
    tmax = 0.0
    time, bodyarray, number_bodies = io_oberon.io_nbody.read_nbody_datafile(input_file, tmax)    
    
    # First set up the figure, the axis, and the plot element we want to animate            

    xmin = 1.0e30
    xmax = -1.0e30
    ymin = xmin
    ymax = xmax
    amax = ymax
    emax = -1.0
    
    linemax = -1
    colours = []
    sizes = []
    
    for i in range(number_bodies):
        
        # Make arrays to define plot sizes and colours
        # Colour is a function of object mass
        # Colour tables for plots (picked from xkcd graphic)
        # Pick object colour
        
        AUtoEarthRadii = 234814.000942
        
        colours.append(pick_circle_colour(bodyarray[i].radius*AUtoEarthRadii))
        sizes.append(100.0/(1.0-np.log10(bodyarray[i].radius)))
                
        # Test for minimum and maximum values of x,y,a and e
        xtry = np.amin(bodyarray[i].x)
        if(xtry<xmin): xmin = xtry
        xtry = np.amax(bodyarray[i].x)
        if(xtry>xmax): xmax = xtry

        ytry = np.amin(bodyarray[i].y)
        if(ytry<ymin): ymin = ytry
        ytry = np.amax(bodyarray[i].y)
        if(ytry>ymax): ymax = ytry
        
        atry = np.amax(bodyarray[i].a)
        if(atry>amax): amax = atry
        
        etry = np.amax(bodyarray[i].e)
        if(etry>emax): emax = etry
        
        nlines = len(bodyarray[i].x)
        
        if(nlines>linemax): linemax = nlines
        
    
    xmin = xmin*0.99
    xmax = xmax*1.01
    ymin = ymin*0.99
    ymax = ymax*1.01
    amax = amax*1.01

    # Create a multiple plot figure
        
    plt.ion()    
    fig2 = plt.figure()

    ax_pos = fig2.add_subplot(211)
    ax_orb = fig2.add_subplot(212)    
    
    npoints = 100
    
    for j in range(linemax):
        
        # Start with x-y plot
        
        ax_pos.set_xlim(xmin,xmax)
        ax_pos.set_ylim(ymin,ymax)
        ax_pos.set_xlabel('x (AU)')
        ax_pos.set_ylabel('y (AU)')
                
        for i in range(number_bodies):
            if(i!=0):
                xorb, yorb, zorb = calc_orbit_track(bodyarray[i].a[j], bodyarray[i].e[j], bodyarray[i].i[j],bodyarray[i].longascend[j],bodyarray[i].argper[j],npoints)
                xorb = xorb-bodyarray[0].x[j]
                yorb = yorb-bodyarray[0].y[j]
                
                if(np.amin(xorb) < xmin): xmin= np.amin(xorb)
                if(np.amax(xorb) > xmax): xmax= np.amax(xorb)
                
                if(np.amin(yorb) < ymin): ymin= np.amin(yorb)
                if(np.amax(yorb) > ymax): ymax= np.amax(yorb)
                
                
                ax_pos.plot(xorb,yorb, color = colours[i])
                
            ax_pos.scatter(bodyarray[i].x[j],bodyarray[i].y[j],facecolor=colours[i], s = sizes[i])
            ax_pos.text(0.6*xmax, 1.1*ymax, 'Time ='+str(time[j])+' yr')
            
        # Now do a-e plot
        ax_orb.set_xlim(0.0,amax)
        ax_orb.set_ylim(0.0,emax)
        ax_orb.set_xlabel('a (AU)')
        ax_orb.set_ylabel('e')
        ax_orb.set_xscale('log')
        #ax_orb.set_yscale('log')
        for i in range(1,number_bodies):
            ax_orb.scatter(bodyarray[i].a[j],bodyarray[i].e[j], facecolor=colours[i], s = sizes[i])                    
                
        sleep(0.01)
        plt.draw()
        ax_pos.clear()
        ax_orb.clear()
        
        
def pick_circle_colour(rad):
    '''Selects the plotting colour for the object depending on its radius'''    

    print rad
    # Sub Earths (lightblue)
    subearth = (29.0/256.0,166.0/256.0,97.0/256.0)
    # Earths
    earth = (20.0/256.0,107.0/256.0,135.0/256.0)
    # Super Earths
    superearth = (135.0/256.0,99/256.0,21/256.0)
    # Neptunes
    neptunes = (84.0/256.0,53.0/256.0,16.0/256.0)
    # Jupiters
    jupiters = (135.0/256.0,40.0/256.0,21.0/256.0)
    # Stars are yellow
    stars = (255.0/256.0, 255.0/256.0, 0.0)
    
    colors = jupiters
                
    # Sub Earth
    if rad < 0.8:
        colors = subearth        
    elif rad >= 0.8 and rad <1.25:
        colors = earth
    # Super Earths
    elif rad >=1.25 and rad < 2.6:
        colors = superearth
    # Neptunes
    elif rad >=2.6 and rad < 6.0:
        colors = neptunes
    # Jupiters
    elif rad >=6.0 and rad < 20.0:
        colors = jupiters
    elif rad>20.0:
        colors = stars 
        
    print colors
    return colors



def calc_orbit_track(semimaj,ecc,inc,longascend,argper, npoints):
    '''Given an input body's orbital parameters, 
    calculates x and y coordinates for
    its orbit over N points'''
        
    if(ecc <1.0):
        nu = np.linspace(0,2.0*np.pi, num=npoints)
    else:
        nu = np.linspace(-np.arccos(1.0/ecc), np.arccos(1.0/ecc), num=npoints)                
    
    # EDIT: get semilat directly from 
    
    semilat = 1.0
    
    if(ecc<1.0):
        semilat = semimaj*(1.0-ecc*ecc)
    elif(ecc>1.0):
        semilat = -semimaj*(1.0-ecc*ecc)    
    
    r = semilat/(1.0+ecc*np.cos(nu))
    
    x = r*(np.cos(longascend)*np.cos(argper+nu) - np.sin(longascend)*np.sin(argper+nu)*np.cos(inc))
    y = r*(np.sin(longascend)*np.cos(argper+nu) - np.cos(longascend)*np.sin(argper+nu)*np.cos(inc))
    z = r*(np.sin(argper+nu)* np.sin(inc))
    
    return x,y,z    
    
        
        
    
    
    
    
if __name__ == '__main__':
    show_orbit()
