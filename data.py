#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis

""" Module data.py

This module define the classes SeismicData() to handle seismic data set.
These datasets define the geographic repartition of raypath in the inner core.

functions:
    read_from_file: to read a file with seismic data

classes:
    SeismicData: base class
    SeismicFromFile: data obtained from a file (real data)
    PerfectSamplingEquator
    PerfectSamplingEquatorRadial
    RandomData: random data, well partitioned on the horizontal, and between 15 and 106km
"""
from __future__ import division
from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt  # for figures
# from mpl_toolkits.basemap import Basemap  # to render maps
import pandas as pd

# personal routines
import positions
# import geodynamic
import plot_data


def read_from_file(filename, names=["station",
                                    "PKIKP-PKiKP travel time residual",
                                    "zeta",
                                    "epicentral distance",
                                    "station lat",
                                    "station lon",
                                    "event lat",
                                    "event lon",
                                    "event depth",
                                    "in lat",
                                    "in lon",
                                    "out lat",
                                    "out lon",
                                    "turn lat",
                                    "turn lon",
                                    "turn depth",
                                    "inner core travel time",
                                    "PKIKP/PKiKP amplitude ratio"],
                   slices="all"):
    """ read seismic data repartition

    input parameters:
    - filename: name of the data file
    - names: names of the columns for the data set
    - slices: names of columns for the output.
    output:
    - data : pandas DataFrame with all the datas.
    Columns name are indicated by the variable "names".
    """
    data = pd.read_table(filename, sep=' ', names=names, skiprows=0)
    if slices != "all":
        data = data[slices]
    return data


class SeismicData(object):
    """ Class for seismic data """

    def __init__(self):
        self.data_points = []
        self.size = None
        self.proxy = 0.
        self.name = None

    def __getitem__(self, key):
        return self.data_points[key]

    def extract_xyz(self, type_of_point):
        """Extract the cartesian coordinates of the points in the data set"""
        x, y, z = np.empty([self.size, 1]), np.empty(
            [self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            point = getattr(ray, type_of_point)
            x[i] = point.x
            y[i] = point.y
            z[i] = point.z
        return x, y, z

    def extract_rtp(self, type_of_point):
        """Extract the radius, theta (latitute), phi (longitude) of the points"""
        r, theta, phi = np.empty([self.size, 1]), np.empty(
            [self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            point = getattr(ray, type_of_point)
            r[i] = point.r
            theta[i] = point.theta
            phi[i] = point.phi
        return r, theta, phi

    def map_plot(self, geodyn_model=''):
        """ plot data on a map."""
        # user should plot on map in the main code.

        m, fig = plot_data.setting_map()
        colormap = plt.cm.get_cmap('RdYlBu')
        _, theta, phi = self.extract_rtp("bottom_turning_point")
        x, y = m(phi, theta)
        proxy = np.array([self.proxy]).T.astype(float)
        sc = m.scatter(x, y, c=proxy, zorder=10, cmap=colormap)

        # TO DO : make a function to plot great circles correctly!
        #r1, theta1, phi1 = self.extract_in() #use extract_rtp()
        #r2, theta2, phi2 = self.extract_out()
        # for i, t in enumerate(theta1):
        #    z, w = m.gcpoints(phi1[i], theta1[i], phi2[i], theta2[i], 200)#
        #    m.plot(z, w, zorder=5, c="black")
        #    m.drawgreatcircle(phi1[i], theta1[i], phi2[i], theta2[i], zorder=5, c="black")
        title = "Dataset: {},\n geodynamic model: {}".format(
            self.name, geodyn_model)
        plt.title(title)
        plt.colorbar(sc)
        # plt.show()

    def phi_plot(self, geodyn_model=''):
        """ Plot proxy as function of longitude """
        # user should use pyplot.plot functions in the main code
        fig, ax = plt.subplots()
        _, _, phi = self.extract_rtp("bottom_turning_point")
        ax.plot(phi, self.proxy, '.')
        title = "Dataset: {},\n geodynamic model: {}".format(
            self.name, geodyn_model)
        plt.title(title)
        plt.xlabel("longitude of bottom turning point")
        plt.ylabel("proxy")
        # plt.show()

    def distance_plot(self, geodyn_model='', point=positions.SeismoPoint(1., 0., 0.)):
        """ Plot proxy as function of the angular distance with point G """
        # user should use pyplot.plot functions in the main code
        fig, ax = plt.subplots()
        _, theta, phi = self.extract_rtp("bottom_turning_point")
        theta1, phi1 = point.theta, point.phi
        distance = positions.angular_distance_to_point(
            theta, phi, theta1, phi1)
        ax.plot(distance, self.proxy, '.')
        title = "Dataset: {},\n geodynamic model: {}".format(
            self.name, geodyn_model)
        plt.title(title)
        plt.xlabel(
            "Angular distance between turning point and ({} {})".format(theta1, phi1))
        plt.ylabel("proxy")
        # plt.show()


class SeismicFromFile(SeismicData):
    """Seismic data set from file."""

    def __init__(self, filename="results.dat", RICB=1221.):
        SeismicData.__init__(self)
        self.name = "Data set from Lauren's file"
        # seismic data set (from Lauren's file)
        self.filename = filename
        self.slices = ["PKIKP-PKiKP travel time residual", "turn lat",
                       "turn lon", "turn depth", "in lat", "in lon", "out lat", "out lon"]
        self.data = read_from_file(filename, slices=self.slices)
        self.size = self.data.shape[0]
        self.data_points = []
        for _, row in self.data.iterrows():
            ray = positions.Raypath()
            ray.add_b_t_point(positions.SeismoPoint(
                1. - row["turn depth"] / RICB, row["turn lat"], row["turn lon"]))
            in_point = positions.SeismoPoint(1., row["in lat"], row["in lon"])
            out_point = positions.SeismoPoint(
                1., row["out lat"], row["out lon"])
            ray.add_in_out(in_point, out_point)
            ray.residual = row["PKIKP-PKiKP travel time residual"]
            self.data_points = np.append(self.data_points, ray)

    def real_residual(self):
        """Extract the values of residuals from the file"""
        value = []
        for ray in self.data_points:
            value = np.append(value, ray.residual)
        return value


class PerfectSamplingEquator(SeismicData):

    def __init__(self, N, rICB=1.):
        SeismicData.__init__(self)
        self.rICB = rICB
        self.N = N
        self.name = "Perfect sampling in the equatorial plane"
        for x in np.linspace(-self.rICB, self.rICB, N):
            for y in np.linspace(-self.rICB, self.rICB, N):
                ray = positions.Raypath()
                ray.add_b_t_point(positions.CartesianPoint(x, y, 0.))
                if ray.bottom_turning_point.r <= self.rICB:
                    self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)

    def plot_c_vec(self, modelgeodyn, proxy=1, cm=plt.get_cmap('summer'), nameproxy=""):
        """ Plot contourf of the proxy + streamlines of the flow.

        Args:
            modelgeodyn: a geodyn.Model instance
            proxy: the values to be plot are either defined as self.proxy,
            given as proxy in the function, or set to 1 if not given.
        """

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        x1 = np.linspace(-self.rICB, self.rICB, self.N)
        y1 = np.linspace(-self.rICB, self.rICB, self.N)
        X, Y = np.meshgrid(x1, y1)
        Z = -1. * np.ones_like(X)
        x, y, _ = self.extract_xyz("bottom_turning_point")
        for it, pro in enumerate(proxy):
            ix = [i for i, j in enumerate(x1) if j == x[it]]
            iy = [i for i, j in enumerate(y1) if j == y[it]]
            Z[ix, iy] = pro
        mask_Z = Z == -1
        Z = np.ma.array(Z, mask=mask_Z)
        sc = ax.contourf(Y, X, Z, 10, cmap=cm)
        #sc2 = ax.contour(Y, X, Z, 10, colors='w')

        Vx, Vy = np.empty((self.N, self.N)), np.empty((self.N, self.N))
        for ix, _ in enumerate(x1):
            for iy, _ in enumerate(y1):
                velocity = modelgeodyn.velocity(
                    modelgeodyn.tau_ic, [X[ix, iy], Y[ix, iy], 0.])
                Vx[ix, iy] = velocity[0]
                Vy[ix, iy] = velocity[1]
        Vx = np.ma.array(Vx, mask=mask_Z)
        Vy = np.ma.array(Vy, mask=mask_Z)
        #ax.quiver(X, Y, Vx, Vy)
        ax.streamplot(X, Y, Vx, Vy, color='black',
                      arrowstyle='->', density=0.5)
        theta = np.linspace(0., 2 * np.pi, 1000)
        ax.plot(np.sin(theta), np.cos(theta), 'k', lw=3)
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])

        cbar = plt.colorbar(sc)
        cbar.set_label(nameproxy)
        title = "Geodynamical model: {}".format(modelgeodyn.name)
        plt.title(title)
        plt.axis("off")
        # plt.show()


class RandomData(SeismicData):
    """ Random repartition of point, between depth 106 and 15km"""

    def __init__(self, N, rICB=1.):
        SeismicData.__init__(self)
        self.rICB = rICB
        self.N = N
        self.name = "Random repartition of data, between 0 and 100km depth"
        self.random_method = "uniform"
        self.depth = [15. / 1221., 106. / 1221.]

        for _ in range(N):
            ray = positions.Raypath()
            ray.add_b_t_point(positions.RandomPoint(
                self.random_method, self.depth, rICB))
            self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)


class PerfectSamplingEquatorRadial(SeismicData):
    """ Points in the equatorial cross section, along radius"""

    def __init__(self, Nr, Ntheta, rICB=1.):
        SeismicData.__init__(self)
        self.rICB = rICB
        self.name = "Perfect sampling in the equatorial plane"
        theta = 0.  # latitude
        for phi in np.linspace(0., 360., Ntheta):
            for r in np.linspace(0.1 * self.rICB, self.rICB * 0.99, Nr):
                ray = positions.Raypath()
                ray.add_b_t_point(positions.SeismoPoint(r, theta, phi))
                self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)

    def radius_plot(self, geodyn_model):
        """ Plot proxy as function of radius (to check growth rate) """
        fig, ax = plt.subplots()
        r, _, phi = self.extract_rtp("bottom_turning_point")
        ax.scatter(r, self.proxy, c=phi, cmap="flag")
        if geodyn_model.proxy_type == "age":
            ax.plot(r, 1.e-6 * geodyn_model.time_unit * (1 - r**2), 'x')
        elif geodyn_model.proxy_type == "growth rate":
            ax.plot(r, 0.5 * geodyn_model.length_unit /
                    geodyn_model.time_unit / r)
        title = "Dataset: {},\n geodynamic model: {}".format(
            self.name, geodyn_model)
        plt.title(title)
        plt.xlabel("radius of point")
        plt.ylabel("proxy")
