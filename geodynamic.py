#!/usr/local/bin/python
# Time-stamp: <2016-02-23 09:44:14 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

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


class Point_evolution():
    """ Geodynamical modelling of a flow. """

    def __init__(self, Point, time=0.):
        """ Point is a positions.Point """
        self.initial_position = Point
        self.old = []
        self.new = Point
        self.time = time
        self.inside = True
        self.history = []

    def advection(self):
        """ advection of the point. Need to be implemented in the derived classes """
        raise NotImplementedError("need to implement advection() in derived class!")

    def check_inside(self):
        pass

    def translation(self, velocity, dt):
        """ translation of the point.

        velocity: tupple with 3 components (cartesian coordinates. x,y,z)
        """
        self.new.x = self.old.x + velocity*dt

class Translation(Point_evolution):

    ### the __init__ method is directly obtained from the class Point_evolution, so no need to redefined it in a derived class
    ## def __init__(self, Point, time)
    ##     Point_evolution.__init__(self, Point, time)

    def advection(self, velocity):
        """ advection of the point. """
        self.old = self.new
        translation(self, velocity, dt)

    def analytical_singlepoint(self, velocity, direction=positions.CartesianPoint(1,0,0)):
        self.exact_solution = exact_translation(self.initial_position, velocity, direction)
        assert(self.exact_solution>0)

    def analytical_raypath(self, velocity, direction):
        pass


if __name__ == '__main__':

    
    
    A = positions.CartesianPoint(1200,0,0)
    Point_evolution(A, 1)
    B = Translation(A,1)
    B.analytical(2000, 1)
    
