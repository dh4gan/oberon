'''
Created on Apr 16, 2013

@author: Tripod
'''

import numpy as np
import csv as c     
from matplotlib import pyplot as py 

class body(object):
    
    

    
    '''
    classdocs
    '''


    def __init__(self, bodytype, mass, radius, x, y, z, vx, vy, vz):
        '''
        Constructor
        '''
        self.bodytype = bodytype
        self.mass = np.float(mass)
        self.radius = np.float(radius)
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.vx = np.array(vx)
        self.vy = np.array(vy)
        self.vz = np.array(vz)
        
    def update_body(self, x, y, z, vx, vy, vz):
        self.x = np.append(self.x, np.float(x))
        self.y = np.append(self.y, np.float(y))
        self.z = np.append(self.z, np.float(z))
        self.vx = np.append(self.vx, np.float(vx))
        self.vy = np.append(self.vy, np.float(vy))
        self.vz = np.append(self.vz, np.float(vz))

        