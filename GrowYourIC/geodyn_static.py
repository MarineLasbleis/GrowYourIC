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
        self.tau_ic = 0.

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
    """ Simple radial symmetry (no flow)

    An additional function can be added and given as a variable (fonction) to define the radial dependency.
    See the class Innermost_IC for an example.
    """

    def __init__(self, fonction=None):
        self.name = "Radial symmetry"
        self.tau_ic = 0.

        if fonction == None:   #TODO maybe should be something as "check if self.radial_dependency is defined, and if not, then defines it"?
            def fonction(r):
                return r
        self.radial_dependency = fonction #has to be a function with 1 argument (radius)

    def proxy_singlepoint(self, point, proxy_type):
        """  """
        proxy = {}  # empty dict
        proxy["radius"] = self.radial_dependency(point.r)
        #np.sign(np.sin((point.phi + self.anglephi) * np.pi / 180.))
        return proxy

    def velocity(self, time, point):
        """ has to be defined, but set to 0 """
        return [0., 0., 0.]

    def radius_ic(self, t):
        return self.rICB

    def verification(self):
        pass


class Innermost_IC(Radial_sym):
    """ """

    def __init__(self, radius_IIC):
        self.radius_IIC = radius_IIC 

        def fonction(r):
            if r>self.radius_IIC:
                answer = 1.
            else:
                answer = 0.
            return answer

        Radial_sym.__init__(self, fonction)
        self.name = "Innermost inner core"
