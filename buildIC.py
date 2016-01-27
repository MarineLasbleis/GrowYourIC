#!/usr/local/bin/python
# Time-stamp: <2016-01-27 16:40:28 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Subproject : Building IC structure
# Author : Marine Lasbleis


import numpy as np
import matplotlib.pyplot as plt

import seismo

def advection(old_position, scheme="none"):
    """ """

    new_position = old_position
    
    return new_position


def growth(old_position, rIC_t, dr, last_number, N):
    """ """
    
    radius_newpoints = rIC_t+dr
    print radius_newpoints
    phi = np.linspace(0, 360, N)
    theta = np.zeros((N))
    add_on = np.concatenate((radius_newpoints*np.ones((N)), theta, phi,
                             np.arange(last_number, last_number+N))).reshape(4, N)
    print np.shape(add_on), np.shape(old_position)
    new_position = np.concatenate((old_position, add_on), axis=1)

    
    return new_position



def insideIC(position, rIC_t):
    """ """
    inside = 1 #1 if true, 0 is false

    return inside


if __name__ == '__main__':

    N = 10 #nb of timesteps
    tau_ic = 1.
    t = np.linspace(0.1,1.,N)

    innercore = [[0], [0], [0]] # each line is : r, theta (latitude), phi, proxy

    # initialisation 

    r = np.linspace(0.1, 0.8, 2)
    #theta = np.linspace(-90, 90, 20)
    phi = np.linspace(0, 360, 20)
    print np.shape(phi)
    for i, j in enumerate(r):
        add_on = np.concatenate((j*np.ones((20)), np.zeros((20)), phi)).reshape(3, 20)
        print np.shape(add_on), np.shape(innercore)
        innercore = np.concatenate((innercore, add_on), axis=1)

    dummy, N_points = np.shape(innercore)
    innercore = np.concatenate((innercore, np.arange(N_points).reshape(1,N_points)), axis=0)

    innercore = advection(innercore)

    innercore = growth(innercore, 0.8, 0.1, N_points, 20)    
        
    
    print innercore
    print np.shape(innercore)
    print np.diff(innercore[3,:])
