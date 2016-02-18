#!/usr/local/bin/python
# Time-stamp: <2016-02-18 16:19:21 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Subproject : Building IC structure
# Author : Marine Lasbleis


import numpy as np
import matplotlib.pyplot as plt

import seismo

def advection(old_position, scheme=None):
    """ Advection of the point old_position """

    if scheme==None:
        new_position = old_position

    else:
        pass
        
    return new_position


def growth(old_position, rIC_t, t, N):
    """ Grow your IC by adding points on the layer at radius rIC_t+dr

    N : number of points in the layer
    """
    
    radius_newpoints = rIC_t
    phi = np.linspace(0, 360, N)
    theta = np.zeros((N))
    add_on = np.concatenate((radius_newpoints*np.ones((N)), theta, phi,
                             t*np.ones(N))).reshape(4, N)
    print np.shape(old_position), np.shape(add_on)
    new_position = np.concatenate((old_position, add_on), axis=1)
    return new_position



def insideIC(position, rIC_t):
    """ """
    inside = 1 #1 if true, 0 is false
    
    return inside

def density_points(N0):
    """ density of points at the surface for r=1"""
    return N0/(4.*np.pi)

def numberpoints_surface(r, density):
    """ return the number of points needed to cover the surface at radius 'r' with a density 'density'

    r: radius (radius = 1 correspond to IC now)
    density: density of points at the surface (corresponding to r=1)
    """
    return int(4.*np.pi*r**2.*density)


if __name__ == '__main__':

    N = 10 #nb of timesteps
    tau_ic = 1.
    t = np.linspace(0.1,1.,N)

    innercore = [[0], [0], [0], [0]] # each line is : r, theta (latitude), phi, proxy

    # initialisation 

    r = np.linspace(0.1, 1., 20)
    #theta = np.linspace(-90, 90, 20)
    ## phi = np.linspace(0, 360, 20)
    ## print np.shape(phi)
    ## for i, j in enumerate(r):
    ##     add_on = np.concatenate((j*np.ones((20)), np.zeros((20)), phi)).reshape(3, 20)
    ##     print np.shape(add_on), np.shape(innercore)
    ##     innercore = np.concatenate((innercore, add_on), axis=1)

    ## dummy, N_points = np.shape(innercore)
    ## innercore = np.concatenate((innercore, 0.4*np.ones(N_points).reshape(1,N_points)), axis=0)
    #    innercore = np.concatenate((innercore, np.arange(N_points).reshape(1,N_points)), axis=0)

    #innercore = advection(innercore)


    density = density_points(100)
    for i, radius in enumerate(r):
        innercore = growth(innercore, radius, radius, numberpoints_surface(radius, density))
        print numberpoints_surface(radius, density)
    #    innercore = growth(innercore, 0.6, 0.1, 0.6, numberpoints_surface(0.6, density))
    #    innercore = growth(innercore, 0.7, 0.1, 0.7, numberpoints_surface(0.7, density))  
        
    fig, ax = plt.subplots(1)
    x, y, z = seismo.from_seismo_to_cartesian(innercore[0,:], innercore[1,:], innercore[2,:])
    ax.scatter(x, y, innercore[3,:]*100)


    fig2, ax2 = plt.subplots()
    ax2.plot(innercore[0,:], innercore[3,:], ".")
    
    print np.shape(innercore)
    
    plt.show()
