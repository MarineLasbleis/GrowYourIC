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
import positions
import intersection
import geodyn


class Hemispheres(geodyn.Model):
    """ Static hemispheres: 

        proxy is just defines as -1 in the western hemisphere and +1 in the eastern one."""

    def __init__(self):
        self.name = "Static hemispheres"

    def proxy_singlepoint(self, point, proxy_type):
        """ -1 in western hemisphere, +1 in the eastern hemisphere"""
        proxy = {}  # empty dict
        angle = -30
        proxy["age"] = np.sign(np.sin((point.phi + angle) * np.pi / 180.))
        return proxy

    def velocity(self, time, point):
        """ has to be defined, but set to 0 """
        return [0., 0., 0.]

    def radius_ic(self, t):
        return self.rICB
