#!/usr/local/bin/python
# Time-stamp: <2016-02-22 22:04:34 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

#personal routines
import positions


RICB = 1221.


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
        pass #TODO : to write !!!


class Translation(Point_evolution):

    ### the __init__ method is directly obtained from the class Point_evolution, so no need to redefined it in a derived class
    ## def __init__(self, Point, time)
    ##     Point_evolution.__init__(self, Point, time)

    def advection(self, velocity):
        """ advection of the point. """
        translation(self, velocity, dt)

    def analytical(self, velocity, time):
        self.exact_solution = (-np.sqrt(RICB**2-self.initial_position.y**2-self.initial_position.z**2)-self.initial_position.x) / velocity +time
        assert(self.exact_solution>0)




if __name__ == '__main__':

    
    
    A = positions.Point(1,0,0,"cartesian")
    Point_evolution(A, 1)
    B = Translation(A,1)
    B.analytical(10,1)
    
