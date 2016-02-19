#!/usr/local/bin/python
# Time-stamp: <2016-02-19 16:19:23 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis


import numpy as np

import sys # .float_info import epsilon # to use assert on floating point equivalence


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
        assert(abs(np.sqrt(self.x**2+self.y**2+self.z**2)-self.r)< 2.*sys.float_info.epsilon)

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


    def straigth_trajectory(self, Point1, Point2, N):
        """ Trajectory is a straigth line between Point1 and Point2, with N points.

        Point1, Point2: Position()
        N: integer (number of points on the trajectory)
        
        Use the cartesian coordinates of both points.
        """
        _Points = []
        _vector = [Point2.x-Point1.x, Point2.y-Point1.y, Point2.z-Point1.z]
        _length = np.sqrt(_vector[0]**2+_vector[1]**2+_vector[2]**2)
        for dx in np.linspace(0, 1, N):
            _Points.append(Position(Point1.x+_vector[0]*dx, Point1.y+_vector[1]*dx, Point1.z+_vector[2]*dx, "cartesian"))
        return _Points, _length
        
    def straigth_in_out(self, N):
        """ Trajectory is a straigth line between in and out points, with N points. """
        if not (self.in_point == None or self.out_point == None):
            self.points = []
            self.points, self.length = self.straigth_trajectory(self.in_point, self.out_point, N)
        else:
            raise Exception("in and out points have not been defined!")
        
    def straigth_in_out_bt(self, N):
        """ Trajectory is a straigth line between in and out points, with 2N-1 points. """
        if not (self.in_point == None or self.out_point == None or self.bottom_turning_point == None):
            points1, length1 = self.straigth_trajectory(self.in_point, self.bottom_turning_point, N)
            points2, length2 = self.straigth_trajectory(self.bottom_turning_point, self.out_point, N)
            self.points = []
            self.length = length1+length2
            self.points = points1 + points2[1:]
            
        else:
            raise Exception("in, out or bottom turning points have not been defined!")


        
if __name__ == '__main__':


    ## position1 = Position(1 , 0, 0, "seismo")
    ## print position1.x, position1.y, position1.z
    ## print position1.r, position1.theta, position1.phi

    ## position1.random_point()
    ## print position1.x, position1.y, position1.z
    ## print position1.r, position1.theta, position1.phi


    ## trajectory = Raypath("BT-point", 1, 0, 0)
    ## print trajectory.bottom_turning_point.r
    ## trajectory.b_t_point()


    trajectory = Raypath("in-out", 2, 0, 0, 1, 0, 90.)
    print "seismo, in:", trajectory.in_point.r, trajectory.in_point.theta, trajectory.in_point.phi
    print "seismo, out:",trajectory.out_point.r, trajectory.out_point.theta, trajectory.out_point.phi
    print "cart, in:",trajectory.in_point.x, trajectory.in_point.y, trajectory.in_point.z
    print "cart, out:", trajectory.out_point.x, trajectory.out_point.y, trajectory.out_point.z
    trajectory.straigth_in_out(10)
    for i in range(10):
        print "seismo", trajectory.points[i].r, trajectory.points[i].theta, trajectory.points[i].phi
    for i in range(10):
        print "cart", trajectory.points[i].x, trajectory.points[i].y, trajectory.points[i].z

    print "==="
    trajectory.b_t_point(1,0,0, "seismo")
    trajectory.straigth_in_out_bt(10)
    for i in range(len(trajectory.points)):
        print "seismo", trajectory.points[i].r, trajectory.points[i].theta, trajectory.points[i].phi
    for i in range(len(trajectory.points)):
        print "cart", trajectory.points[i].x, trajectory.points[i].y, trajectory.points[i].z
