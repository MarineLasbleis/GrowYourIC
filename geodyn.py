#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis

""" geodyn.py

Define the base class for geodynamical models

CLass:
    Model
        def set_parameters
        def verification
        def proxy_singlepoint: NotImplemented
        def define_units
Functions:
    evaluate_proxy: evaluate proxy over all points of the data set using the geodyn model.
    average_proxy: method to average proxy over raypath
"""

from __future__ import division
from __future__ import absolute_import


import numpy as np


class Model(object):
    """Base class for all the geodynamical models"""


    def set_parameters(self, dict_param):
        """ add any parameter of the form {'param':value} as self.param = value """
        for k, v in dict_param.items():
            setattr(self, k, v)

    def verification(self):
        """Verify if everything required is here.

        Verify if the geodynamical model verify some very simple assumptions
        (ex: non zero radius, translation velocity fast enough if case without growth, etc.)
        this method has to be implemented in derived class.
        By default, no verification is done."""
        raise NotImplementedError(
            "need to implement verification() in derived class!")

    def proxy_singlepoint(self, point):
        """ evaluate the proxy on a single positions.Point instance."""
        raise NotImplementedError(
            "need to implement proxy_singlepoint() in derived class!")

    def define_units(self):
        """Define the units"""
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
            self.time_unit, self.length_unit = 1e9, 1.221e6
           # if no units have been given, we set them arbitrarly


def evaluate_proxy(dataset, method, proxy="", verbose=True):
    """ evaluate the value of the proxy on all the points of the data set

        dataset : a data.SeismicData object
        method : a geodyn.Model  object
        """
    if proxy=="": proxy = method.proxy_type
    print("===")
    print("== Evaluate value of proxy for all points of the data set ")
    print("= Geodynamic model isi {}".format(method.name))
    print("= Proxy is {}".format(proxy))
    print("= Data set is {}".format(dataset.name))
    print("= Proxy is evaluated for {}".format(dataset.method))
    if dataset.method == "raypath":
        print("=== Raypath is {} number of points".format(dataset.NpointsRaypath))
    print("= Number of points to examine: {}".format(dataset.size))

    method.verification()

    proxy = np.empty_like(dataset.data_points)
    for i, ray in enumerate(dataset.data_points):
        if i % 100 == 0 and verbose:
            print("Computing Ray number {}".format(i))
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
            number_points = dataset.NpointsRaypath
            dataset.data_points[i].straigth_in_out(number_points + 2)
            raypath = ray.points
            proxy[i] = average_proxy(raypath, method)
    print("===")
    return np.array(proxy).astype(float)


def average_proxy(ray, method):
    r""" method to average proxy over the raypath.

    Simple method is direct average of the proxy: $\sum proxy(r) / \sum dr$.
    Other methods could be: $1/(\sum 1 / proxy)$ (better for computing \delta t)
    """
    total_proxy = 0.
    if method.evaluation == "inverse":
        for _, point in enumerate(ray):
            _proxy = method.proxy_singlepoint(point)[method.proxy_type]
            total_proxy += 1. / _proxy
        number = len(ray)
        proxy = 1. / total_proxy / float(number)
    else:
        for j, point in enumerate(ray):
            _proxy = method.proxy_singlepoint(point)[method.proxy_type]
            total_proxy += _proxy
        number = len(ray)
        proxy = total_proxy / float(number)

    return proxy
