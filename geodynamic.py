#!/usr/local/bin/python
# Time-stamp: <2016-02-23 09:44:14 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math

#personal routines
import positions


RICB = 1221.


def exact_translation(point, velocity, direction=positions.CartesianPoint(1,0,0)):
    x_0, y_0, z_0 = point.x, point.y, point.z
    mean_direction = np.sqrt(direction.x**2+ direction.y**2+ direction.z**2)
    a, b, c = direction.x/mean_direction, direction.y/mean_direction, direction.z/mean_direction

    solution_1 = x_0*a+y_0*b+z_0*c + np.sqrt((x_0*a+y_0*b+z_0*c)**2-(x_0**2+y_0**2+z_0**2-RICB**2))
    solution_2 = x_0*a+y_0*b+z_0*c - np.sqrt((x_0*a+y_0*b+z_0*c)**2-(x_0**2+y_0**2+z_0**2-RICB**    2))
    ## TO DO : verify that we can remove solution_2 ? Is solution_1 always max?
    return max(solution_1, solution_2)


def numberpoints_surface(r, density):
    """ return the number of points needed to cover the surface at radius 'r' with a density 'density'
    
    r: radius (radius = 1 correspond to IC now)
    density: density of points at the surface (corresponding to r=1)
    """
    return int(4.*np.pi*r**2.*density)


class Model():
    """ Geodynamical modelling of a flow. """

    def __init__(self):
        """ Point is a positions.Point """
        self.points = []
        self.old = []
        self.new = [] 
        self.time = None
        self.history = []
        self.radius = None
        self.Npoints = None
        self.density = None

    def initialisation(self, t0, R0, density):
        self.time = t0
        self.radius = R0
        self.density = density
        self.fill_sphere(density) #fill the initial sphere with points. If the initial sphere is too small, only the point in r=0 is set up.
        self.Npoints = len(self.points)

    def run(self, tmax, dt):
        t0 = self.time
        for t in np.arange(t0, tmax, dt):
            self.old = self.points
            self.time = t
            self.advect()
            self.grow()
            self.check_inside()
            self.add_points()
            self.Npoints = len(self.points)
            


    def advect(self):
        """ advection of the point. Need to be implemented in the derived classes """
        raise NotImplementedError("need to implement advection() in derived class!")

    def grow(self):
        """ grow the radius of the core """
        raise NotImplementedError("need to implement grow() in derived class!")

    def check_inside(self):
        """ check if the points are inside or outside the core"""
        # find points outside
        index_outside = [] #store the index of points to remove
        for i, point in enumerate(self.points):
            if point.r > rICB:
                index_outside.append(i)
        # remove them (store?)
        self.points = np.delete(self.points, index_outside)

    def extract_rtp(self):
        r, theta, phi = np.empty([self.Npoints,1]), np.empty([self.Npoints,1]), np.empty([self.Npoints,1])
        for i, point in enumerate(self.points):
            r[i] = point.r
            theta[i] = point.theta
            phi[i] = point.phi
        return r, theta, phi

    def extract_xyz(self):
        x, y, z = np.empty([self.Npoints,1]), np.empty([self.Npoints,1]), np.empty([self.Npoints,1])
        for i, point in enumerate(self.points):
            x[i] = point.x
            y[i] = point.y
            z[i] = point.z
        return x, y, z 

    def add_points(self, density_surface):
        """ add points where there are few of them and on the surface if ic is growing """
        # on the surface
        S = 4. * np.pi * self.radius**2 # size of the surface
        n = int(np.ceil(S*density_surface)) #number of points to add on the surface. 
        for i in range(n):
            u, v = np.random.uniform(0, 1), np.random.uniform(0, 1)
            r = self.radius
            theta = 90. - 180./np.pi*np.arccos(2*v-1) #here is the latitude
            phi = 360. * u
            self.points.append(positions.SeismoPoint(r, theta, phi))
        self.Npoints = len(self.points)


    def fill_sphere(self, density):
        volume = 4./3. * np.pi *self.radius**3. #volume of the sphere
        n = int(math.ceil(density*volume))
        if n == 1:
            self.points.append(positions.SeismoPoint(0.,0.,0.)) #central point
        else:
            # to get a uniform repartition, make it uniform in the (x, y, z) space, but with more points that needed, and then remove all the extra points. Check if at least one point is in the space, and if not add a point at the center. 
            n = int(math.ceil(density*volume*6./np.pi))
            print n
            for i in range(n):
                P = positions.CartesianPoint(self.radius*np.random.uniform(-1., 1.), self.radius*np.random.uniform(-1., 1.), self.radius*np.random.uniform(-1., 1.))
                if P.r <= self.radius:
                    self.points.append(P)
            self.Npoints = len(self.points)
            print self.Npoints
            if self.Npoints <= 1:
                self.points = np.array([positions.SeismoPoint(0.,0.,0.)])
                self.Npoints = 1
                

    #def translation(self, velocity, dt):
    #    """ translation of the point.
    #
    #         velocity: tupple with 3 components (cartesian coordinates. x,y,z)
    #       """
    #       self.new.x = self.old.x + velocity*dt

class Translation(Model):
    ## Not up-to-date
    ### the __init__ method is directly obtained from the class Point_evolution, so no need to redefined it in a derived class
    ## def __init__(self, Point, time)
    ##     Point_evolution.__init__(self, Point, time)

    def advect(self, velocity):
        """ advection of the point. """
        self.old = self.new
        translation(self, velocity, dt)

    def analytical_singlepoint(self, velocity, direction=positions.CartesianPoint(1,0,0)):
        self.exact_solution = exact_translation(self.initial_position, velocity, direction)
        assert(self.exact_solution>0)

    def analytical_raypath(self, velocity, direction):
        pass


if __name__ == '__main__':

    
    
    A = Model()
    print A.Npoints
    A.initialisation(1., 10.,1.) 
    print A.Npoints
    x, y, z = A.extract_xyz()
    plt.plot(np.sqrt(x**2+y**2+z**2), '.')
    plt.show()

    A.add_points(10.)

    print A.Npoints

    x, y, z = A.extract_xyz()
    plt.plot(np.sqrt(x**2+y**2+z**2), '.')
    plt.show()

    r, t, p = A.extract_rtp()
    plt.plot(np.cos(np.pi/2-t*np.pi/180.), '.')
    plt.show()
    plt.plot(p, '.')
    plt.show()
