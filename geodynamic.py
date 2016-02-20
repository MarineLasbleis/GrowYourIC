#!/usr/local/bin/python
# Time-stamp: <2016-02-20 23:18:21 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

#personal routines
import positions



class Point_evolution():
    """ Geodynamical modelling of a flow. """

    def __init__(self, Point, time):
        """ Point is a positions.Point """
        self.initial_position = Points
        self.old = []
        self.new = Points
        self.time = time
        self.inside = True


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
