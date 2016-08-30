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
#import positions
#import intersection


class Model():

    def set_parameters(self, dict_param):
        """ add any parameter of the form {'param':value} as self.param = value """
        for k, v in dict_param.items():
            setattr(self, k, v)

    def verification(self):
        """ verify if the geodynamical model verify some very simple assumptions (ex: non zero radius, translation velocity fast enough if case without growth, etc.) 
        this method has to be implemented in derived class.
        By default, no verification is done."""
        self.units()
        print "No verification has been implemented for this geodynamical model. If needed, please implement them in the class."

    def proxy_singlepoint(self, point):
        """ evaluate the proxy on a single positions.Point instance."""
        raise NotImplementedError(
            "need to implement proxy_singlepoint() in derived class!")

    def define_units(self):
        # parameters have been given in non-dimensional forms.
        if self.units == None:
            # if units have been given separately
            # we need (time_unit, length_unit), (time_unit, velocity_unit) or
            # (length_unit, velocity_unit). For now, let's say we have
            # time_unit and length_unit.
            try:
                self.time_unit
                self.length_unit
            except AttributeError:  # at least one of them does not exist
                # age inner core of 1billion years. age is in years (please be
                # careful)
                self.time_unit = 1.e9
                # radius of inner core today, in km.
                self.length_unit = 1.221e6
        # self.units has to be of the form ("name", "name"), with the two names
        # the useful values to non dimensionalise
        elif self.units == ("tau_ic", "rICB"):
            self.time_unit, self.length_unit = self.tau_ic, self.rICB
        elif self.units == ("vt", "rICB"):
            self.time_unit, self.length_unit = self.rICB / self.vt, self.rICB
        else:
            time_unit, length_unit = 1e9, 1.221e6

           # if no units have been given, we set them arbitrarly


def evaluate_proxy(dataset, method, verbose=True):
    """ evaluate the value of the proxy on all the points of the data set, using the choosen geodynamical method

        dataset : a data.SeismicData object
        method : a geodyn.Model  object
        """
    print "==="
    print "== Evaluate value of proxy for all points of the data set "
    print "= Geodynamic model is", method.name
    print "= Proxy is", method.proxy_type
    print "= Data set is", dataset.name
    print "= Proxy is evaluated for", dataset.method
    if dataset.method == "raypath":
        print "=== Raypath is", dataset.NpointsRaypath, " number of points"
    print "= Number of points to examine: ", dataset.size

    method.verification()

    proxy = np.empty_like(dataset.data_points)
    for i, ray in enumerate(dataset.data_points):
        if i % 100 == 0 and verbose:
            print "Computing Ray number", i
        if dataset.method == "bt_point":
            point = ray.bottom_turning_point
            proxy[i] = method.proxy_singlepoint(point)[method.proxy_type]
        elif dataset.method == "raypath":
            # the raypath is given with constant intervals between points.
            # the average is directly sum()/number of points.
            # ATTENTION: here we compute the average of the proxy
            # but if we want to compute the velocity of Pwaves,
            # we may want to compute the velocity for each point
            # and average over the inverse of the velocities.
            N = dataset.NpointsRaypath
            dataset.data_points[i].straigth_in_out(N + 2)
            raypath = ray.points
            proxy[i] = average_proxy(raypath, method)
    print "==="
    return np.array(proxy).astype(float)


def average_proxy(ray, method):
    """ method to average proxy over the raypath.

    Simple method is direct average of the proxy: \sum proxy(r) / \sum dr.
    Other methods could be: 1/(\sum 1/proxy) (better for computing \delta t)
    """
    total_proxy = 0.
    if method.evaluation == "inverse":
        for j, point in enumerate(ray):
            _proxy = method.proxy_singlepoint(point)[method.proxy_type]
            total_proxy += 1. / _proxy
        N = len(ray)
        proxy = 1. / total_proxy / float(N)
    else:
        for j, point in enumerate(ray):
            _proxy = method.proxy_singlepoint(point)[method.proxy_type]
            total_proxy += _proxy
        N = len(ray)
        proxy = total_proxy / float(N)

    return proxy
