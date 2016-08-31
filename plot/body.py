'''
Created on Apr 16, 2013

@author: Tripod
'''

import numpy as np
 

class body(object):
    
    

    
    '''
    classdocs
    '''


    def __init__(self, bodytype, mass, radius, xi, yi, zi, vxi, vyi, vzi, a,e,i,longascend,argper,meananom):
        '''
        Constructor
        '''
        self.bodytype = bodytype
        self.mass = np.float(mass)
        self.radius = np.float(radius)
        self.x = np.array(float(xi))
        self.y = np.array(float(yi))
        self.z = np.array(float(zi))
        
        self.r = np.sqrt(float(xi)*float(xi) + float(yi)*float(yi) +float(zi)*float(zi))
        self.vx = np.array(float(vxi))
        self.vy = np.array(float(vyi))
        self.vz = np.array(float(vzi))
        self.v = np.sqrt(float(vxi)*float(vxi) + float(vyi)*float(vyi) +float(vzi)*float(vzi))
        
        self.a = np.array(float(a))
        self.e = np.array(float(e))
        self.i = np.array(float(i))
        self.longascend = np.array(float(longascend))
        self.argper = np.array(float(argper))
        self.meananom = np.array(float(meananom))
        
    def update_body(self, xi, yi, zi, vxi, vyi, vzi,a,e,i,longascend,argper,meananom):
        self.x = np.append(self.x, np.float(xi))
        self.y = np.append(self.y, np.float(yi))
        self.z = np.append(self.z, np.float(zi))
        self.r = np.append(self.r, np.sqrt(float(xi)*float(xi) + float(yi)*float(yi) +float(zi)*float(zi)))
        self.vx = np.append(self.vx, np.float(vxi))
        self.vy = np.append(self.vy, np.float(vyi))
        self.vz = np.append(self.vz, np.float(vzi))
        
        
        self.v = np.append(self.v,np.sqrt(float(vxi)*float(vxi) + float(vyi)*float(vyi) +float(vzi)*float(vzi)))
        self.a = np.append(self.a,np.float(a))
        self.e = np.append(self.e,np.float(e))
        self.i = np.append(self.i, np.float(i))
        self.longascend = np.append(self.longascend, np.float(longascend))
        self.argper = np.append(self.argper, np.float(argper))
        self.meananom = np.append(self.meananom, np.float(meananom))

        
