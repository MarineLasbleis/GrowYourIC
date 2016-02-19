#!/usr/local/bin/python
# Time-stamp: <2016-02-19 11:19:50 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis


import numpy as np


def from_seismo_to_cartesian(r, theta, phi):
    """ Calculate the cartesian coordinates from spherical (w/ latitude) coordinates)

    input: 
    r : radius (km)
    theta : latitude (degree)
    phi : longitude (degree)
    output:
    x, y, z
    
    """
    theta = (90-theta) * np.pi/180. # colatitude in rad
    phi = phi * np.pi/180.
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z


def from_cartesian_to_seismo(x, y, z):
    """ Calculate the spherical coordinates (w/ latitude) from cartesian coordinates)

    r, theta, phi = from_cartesian_to_seismo(x, y, z)
    input: x, y, z (same length)
    output:
    r : radius (km)
    theta : latitude (degree)
    phi : longitude (degree)
    (same length as the input)

    """
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arccos(z/r)*180./np.pi # colatitude, in degree
    if (x, y) == (0, 0):
        phi = 0.
    else:
        phi = np.where(y>=0., np.arccos(x/np.sqrt(x**2+y**2)), 2.*np.pi-np.arccos(x/np.sqrt(x**2+y**2)) )
        phi = phi*180./np.pi
        #phi = np.arctan(y/x)*180./np.pi # longitude, in degree
    return r, 90.-theta, phi


class Position:
    """ Position of a point in the Earth.

    can be computed in cartesian coordinates or in "seismological" coordinates.
    Cartesian coordinates: x,y,z (z is the NS, y is the EW and axe x cross the 0 longitude)
    Seismological coordinates: r, theta, phi (theta is the latitude)
    """

    def __init__(self, a, b, c, set_method):
        if set_method == "cartesian":
            self.x, self.y, self.z = float(a), float(b), float(c)
            self.r, self.theta, self.phi = from_cartesian_to_seismo(a, b, c)
        elif set_method == "seismo":
            self.r, self.theta, self.phi = float(a), float(b), float(c)
            self.x, self.y, self.z = from_seismo_to_cartesian(a, b, c)

    def random_point(self, set_method="uniform"):#,type_="turningpoint", seismo="surface"):
        """ Create a random point (not raypath)

        type: type of the distribution. Default is uniform over the sphere of radius self.r
        """

        self.phi = np.random.uniform(-180., 180.)
        self.theta = (np.arccos(2*np.random.uniform(0.,1.)-1)*180./np.pi)-90
        self.x, self.y, self.z = from_seismo_to_cartesian(self.r, self.phi, self.theta)

        #TODO : set other methods of randomisation!



class Raypath:
    """ Raypath inside Inner Core.

    raypath are defined either by:
    - bottom turning point + direction (useful for random trajectories)
    - in and out points (at the surface of IC, useful if coming from real data set)
    """
    
    def __init__(self, set_method, *arg):
        self.points = None
        self.bottom_turning_point = None
        self.in_point = None
        self.out_point = None
        if set_method == "BT-point":
            if len(arg) == 3:
                arg = arg + ("seismo",) #by default, assuming seismo-like coordinates.
            self.bottom_turning_point = Position(arg[0], arg[1], arg[2], arg[3])
        elif set_method == "in-out": #assume seismo-like coordinate
            self.in_point = Position(arg[0], arg[1], arg[2], "seismo")
            self.out_point = Position(arg[3], arg[4], arg[5], "seismo")

    def b_t_point(self, *arg):
        """ Bottom turning point of the trajectory """
        if  self.bottom_turning_point == None:
            assert(len(arg)==4), "coordinates need to be on the form a,b,c,type"
            assert(arg[3]=="seismo" or arg[3]=="cartesian"), 'types of coordinates not well defined!'
            self.bottom_turning_point = Position(arg[0], arg[1], arg[2], arg[3])
        else:
            print "Bottom point of raypath already calculated"

    def straigth_in_out(self, N):
        """ Trajectory is a straigth line between in and out points, with N points.

        Use the cartesian coordinates of both in and out points.
        """
        if not (self.in_point == None or self.out_point == None):
            self.points = []
            vector = [self.out_point.x-self.in_point.x, self.out_point.y-self.in_point.y, self.out_point.z-self.in_point.z]
            self.length = np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
            for dx in np.linspace(0, 1, N):
                self.points.append(Position(self.in_point.x+vector[0]*dx, self.in_point.y+vector[1]*dx, self.in_point.z+vector[2]*dx, "cartesian"))
        else:
            raise Exception("in and out points have not been defined!")
            
    

        
if __name__ == '__main__':


    position1 = Position(1 , 0, 0, "seismo")
    print position1.x, position1.y, position1.z
    print position1.r, position1.theta, position1.phi

    position1.random_point()
    print position1.x, position1.y, position1.z
    print position1.r, position1.theta, position1.phi


    trajectory = Raypath("BT-point", 1, 0, 0)
    print trajectory.bottom_turning_point.r
    trajectory.b_t_point()


    trajectory = Raypath("in-out", 1, 0, 0, 1, 90, 0)
    print trajectory.bottom_turning_point
    trajectory.b_t_point(0,0,0, "seismo")
    print trajectory.bottom_turning_point.x
    trajectory.straigth_in_out(10)
    print trajectory.points
    for i in range(10):
        print trajectory.points[i].r, trajectory.points[i].theta, trajectory.points[i].phi
