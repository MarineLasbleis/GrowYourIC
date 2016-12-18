#!/usr/local/bin/python
# Time-stamp: <2016-02-22 21:21:17 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis

from __future__ import division
from __future__ import absolute_import



import numpy as np

import sys  # .float_info import epsilon # to use assert on floating point equivalence


def from_seismo_to_cartesian(r, theta, phi):
    """ Calculate the cartesian coordinates from spherical (w/ latitude) coordinates)

    input:
    r : radius (km)
    theta : latitude (degree)
    phi : longitude (degree)
    output:
    x, y, z

    """
    theta = (90 - theta) * np.pi / 180.  # colatitude in rad
    phi = phi * np.pi / 180.
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
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
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r) * 180. / np.pi  # colatitude, in degree
    if (x, y) == (0, 0):
        phi = 0.
    else:
        phi = np.where(y >= 0., np.arccos(x / np.sqrt(x**2 + y**2)),
                       2. * np.pi - np.arccos(x / np.sqrt(x**2 + y**2)))
        phi = phi * 180. / np.pi
        # phi = np.arctan(y/x)*180./np.pi # longitude, in degree
    return r, 90. - theta, phi


def angular_distance_to_point(theta1, phi1, theta2, phi2):
    """ angular distance between the point (theta1, phi1) and the point (theta2, phi2) at the surface of the core.

    Args: 
        theta1, theta2 : latitude (degree)
        phi1, phi2: longitude (degree)
    Return phi: angle between the two points (in degree)
    """
    theta1, phi1, theta2, phi2 = theta1 * np.pi / 180., phi1 * \
        np.pi / 180., theta2 * np.pi / 180., phi2 * np.pi / 180.
    return np.arccos(np.sin(theta1) * np.sin(theta2) + np.cos(theta1) * np.cos(theta2) * np.cos(abs(phi1 - phi2))) * 180. / np.pi


class Point():
    """ Position of a point in the Earth.

    can be computed in cartesian coordinates or in "seismological" coordinates.
    Cartesian coordinates: x,y,z (z is the NS, y is the EW and axe x cross the 0 longitude)
    Seismological coordinates: r, theta, phi (theta is the latitude)
    """

    def __init__(self):

        self.x, self.y, self.z, self.r, self.theta, self.phi = None, None, None, None, None, None
       # if set_method == "cartesian":
       #     self.x, self.y, self.z = float(a), float(b), float(c)
       #     self.r, self.theta, self.phi = from_cartesian_to_seismo(a, b, c)
       # elif set_method == "seismo":
       #     self.r, self.theta, self.phi = float(a), float(b), float(c)
       #     self.x, self.y, self.z = from_seismo_to_cartesian(a, b, c)
        # print "(r,t,p)", self.r, self.theta, self.phi
        # print "depsilon, r, r'", abs(np.sqrt(self.x**2+self.y**2+self.z**2)-self.r), self.r, np.sqrt(self.x**2+self.y**2+self.z**2)
       # assert(abs(np.sqrt(self.x**2+self.y**2+self.z**2)-self.r)< 1e4*sys.float_info.epsilon)

    def add_cartesian(self):
        assert(self.r != None)
        assert(self.phi != None)
        assert(self.theta != None)
        self.x, self.y, self.z = from_seismo_to_cartesian(
            self.r, self.theta, self.phi)

    def add_seismo(self):
        assert(self.x != None)
        assert(self.y != None)
        assert(self.z != None)
        self.r, self.theta, self.phi = from_cartesian_to_seismo(
            self.x, self.y, self.z)

    def dimensionless(self, lengthscale):
        self.r = self.r / lengthscale
        self.x, self.y, self.z = self.x / lengthscale, self.y / \
            lengthscale, self.z / lengthscale

    def er(self):
        """ return the cartesian coordinnates of \vec{e}_r.
        """
        try:
            assert(self.r != None)
            assert(self.phi != None)
            assert(self.theta != None)
        except (AttributeError, NameError, AssertionError):
            self.add_seismo()
        phi = self.phi / 180. * np.pi
        theta = (90. - self.theta) * np.pi / 180.
        return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


    def proj_er(self, vector):
        """ projection of a vector on e_r, the radial vector in spherical coordinates.
    
        input: vector (cartesian coordinates)
        output: scalar
        """
        vx, vy, vz = vector[0], vector[1], vector[2] #cartesian coordinates
        phi = self.phi / 180. * np.pi
        theta = (90. - self.theta) * np.pi / 180.
        return np.sin(theta)*np.cos(phi)*vx+ np.sin(theta)*np.sin(phi)*vy+ np.cos(theta)*vz



    # ,type_="turningpoint", seismo="surface"):
    def random_point(self, set_method="uniform", depth=[0., 1.], rICB=1.):
        """ Create a random point (not raypath)

        type: type of the distribution. Default is uniform over the sphere of radius self.r
RICB = 1221.
        """
        r = rICB - np.random.uniform(depth[0], depth[1])
        phi = np.random.uniform(-180., 180.)
        theta = (np.arccos(2 * np.random.uniform(0., 1.) - 1)
                 * 180. / np.pi) - 90
        return r, theta, phi
#           #TODO : set other methods of randomisation!


class SeismoPoint(Point):

    def __init__(self, a, b, c):
        self.r, self.theta, self.phi = a, b, c
        self.add_cartesian()


class CartesianPoint(Point):

    def __init__(self, a, b, c):
        self.x, self.y, self.z = a, b, c
        self.add_seismo()


class RandomPoint(Point):

    def __init__(self, method, depth, rIC=1.):
        self.r, self.theta, self.phi = self.random_point(method, depth, rIC)
        self.add_cartesian()


class Raypath():
    """ Raypath inside Inner Core.

    raypath are defined either by:
    - bottom turning point + direction (useful for random trajectories)
    - in and out points (at the surface of IC, useful if coming from real data set)
    """

    def __init__(self):
        self.points = None
        self.bottom_turning_point = None
        self.direction = None
        self.in_point = None
        self.out_point = None

        #self.bottom_turning_point = Point(arg[0], arg[1], arg[2], arg[3])
        # elif set_method == "in-out": #assume seismo-like coordinate
        #    self.in_point = Point(arg[0], arg[1], arg[2], "seismo")
        #    self.out_point = Point(arg[3], arg[4], arg[5], "seismo")

    def add_b_t_point(self, point):
        """ Bottom turning point of the trajectory """
        assert(self.bottom_turning_point == None)
        self.bottom_turning_point = point

    def add_direction(self, zeta):
        self.direction = zeta

    def add_in_out(self, point_in, point_out):
        self.in_point = point_in
        self.out_point = point_out

    def straigth_trajectory(self, Point1, Point2, N):
        """ Trajectory is a straigth line between Point1 and Point2, with N points.

        Point1, Point2: Point()
        N: integer (number of points on the trajectory)

        Use the cartesian coordinates of both points.
        """
        _Points = []
        _vector = [Point2.x - Point1.x, Point2.y -
                   Point1.y, Point2.z - Point1.z]
        _length = np.sqrt(_vector[0]**2 + _vector[1]**2 + _vector[2]**2)
        for dx in np.linspace(0, 1, N):
            _Points.append(CartesianPoint(
                Point1.x + _vector[0] * dx, Point1.y + _vector[1] * dx, Point1.z + _vector[2] * dx))
        return _Points[1:-1], _length

    def straigth_in_out(self, N):
        """ Trajectory is a straigth line between in and out points, with N points. """
        try:
            self.points = []
            self.points, self.length = self.straigth_trajectory(
                self.in_point, self.out_point, N)
        except(NameError, AttributeError):
            raise Exception("in and out points have not been defined!")

    def straigth_in_out_bt(self, N):
        """ Trajectory is a straigth line between in and out points, with 2N-1 points. """
        if not (self.in_point == None or self.out_point == None or self.bottom_turning_point == None):
            points1, length1 = self.straigth_trajectory(
                self.in_point, self.bottom_turning_point, N)
            points2, length2 = self.straigth_trajectory(
                self.bottom_turning_point, self.out_point, N)
            self.points = []
            self.length = length1 + length2
            self.points = points1 + points2[1:]
        else:
            raise Exception(
                "in, out or bottom turning points have not been defined!")


class Raypath_BT(Raypath):
    """ Raypath defined primarly by the bottom point """

    def __init__(self, point, zeta):
        Raypath.__init__(self)
        self.add_b_t_point(point)
        self.add_direction(zeta)


class Raypath_inout(Raypath):
    """ Raypath defined primarly by the in and out points """

    def __init__(self, point_in, point_out):
        Raypath.__init__(self)
        self.add_in_out(point_in, point_out)


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    theta = 0.
    phi = np.linspace(0, 360, 50)

    v_t = [0.,1,0]

    fig, ax = plt.subplots()

    for p in phi:
        point = SeismoPoint(1, theta, p)
        print(point.r, point.theta, point.phi)
        growth = point.proj_er(v_t)
        ax.plot(p, growth, 'o')

    plt.show()
