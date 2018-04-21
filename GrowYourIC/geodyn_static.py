#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis


import numpy as np
import matplotlib.pyplot as plt  # for figures
from mpl_toolkits.basemap import Basemap  # to render maps
import math
from scipy.integrate import ode
from scipy.optimize import fsolve

# personal routines
from . import positions
from . import intersection
from . import geodyn


class Hemispheres(geodyn.Model):
    """ Static hemispheres: 

        proxy is just defines as -1 in the western hemisphere and +1 in the eastern one."""

    def __init__(self, angletheta=0., anglephi=30.):
        self.name = "Static hemispheres"
        self.anglephi = anglephi
        self.angletheta = angletheta

    def proxy_singlepoint(self, point, proxy_type):
        """ -1 in western hemisphere, +1 in the eastern hemisphere"""
        proxy = {}  # empty dict
        proxy["hemisphere"] = np.sign(np.sin((point.phi + self.anglephi) * np.pi / 180.))
        return proxy

    def velocity(self, time, point):
        """ has to be defined, but set to 0 """
        return [0., 0., 0.]

    def radius_ic(self, t):
        return self.rICB

    def verification(self):
        pass






class Radial_sym(geodyn.Model):
    """ Static hemispheres: 

        proxy is just defines as -1 in the western hemisphere and +1 in the eastern one."""

    def __init__(self, fonction=None):
        self.name = "Radial symmetry"

        if fonction == None:
            def fonction(r):
                return r
        self.function_radius = fonction #has to be a function with 1 argument (radius)

    def proxy_singlepoint(self, point, proxy_type):
        """  """
        proxy = {}  # empty dict
        proxy["radius"] = self.function_radius(point.r)
        #np.sign(np.sin((point.phi + self.anglephi) * np.pi / 180.))
        return proxy

    def velocity(self, time, point):
        """ has to be defined, but set to 0 """
        return [0., 0., 0.]

    def radius_ic(self, t):
        return self.rICB

    def verification(self):
        pass
