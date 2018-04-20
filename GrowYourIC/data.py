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
    PerfectSamplingSurface
"""
from __future__ import division
from __future__ import absolute_import
#TO DO: verify if this is necessary? check with Python 2?



import numpy as np
import matplotlib.pyplot as plt  # for figures
# from mpl_toolkits.basemap import Basemap  # to render maps
import pandas as pd

# personal routines
from . import positions
# import geodynamic
from . import plot_data


# Read from files with pandas: 
## example: pd.read_table(self.filename, sep='\s+', names=self.slices, skiprows=10) 
## example: pd.read_table(self.filename, sep='\s+', header=None)[nb_slices]





class SeismicData(object):
    """ Class for seismic data """

    def __init__(self):
        self.data_points = []
        self.size = None
        self.proxy = 0.
        self.name = None
        self.shortname = None

    def __getitem__(self, key):
        return self.data_points[key]

    def extract_xyz(self, type_of_point="bottom_turning_point"):
        """Extract the cartesian coordinates of the points in the data set"""
        x, y, z = np.empty([self.size, 1]), np.empty(
            [self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            point = getattr(ray, type_of_point)
            x[i] = point.x
            y[i] = point.y
            z[i] = point.z
        return x, y, z

    def extract_rtp(self, type_of_point="bottom_turning_point"):
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


class SeismicFromFile(SeismicData):
    """Seismic data set from file."""

    def __init__(self, filename="WD11.dat", RICB=1221., name="Data set from Waszek and Deuss 2011", shortname="WD11", N="all"):
        SeismicData.__init__(self)
        self.filename = filename
        self.rICB = RICB
        if N=="all": 
            self.limited_nber_points = [False, 0]
        else: 
            self.limited_nber_points = [True, N]
        self.isitknowndataset()

    def isitknowndataset(self, verbose=True):
        """ Check if the data set is already known. If not, explain how to add one. 

        Required variables to specify: 
        self.name and self.shortname : names to be printed on figures and filenames (text)
        self.data_points : all the raypaths (numpy array)
        self.size : total size of the data set (int, number of points)

        """
        if self.filename[-8:] == "WD11.dat":
            self.name = "Data set from Waszek and Deuss 2011"
            self.shortname = "WD11"
            self.WD11()
            #if verbose: 
            print("Waszek and Deuss 2011 successfully loaded. {} trajectories.".format(self.size))
        elif self.filename[-24:] == "DF_sample_ksi_sorted.dat":
            self.name = "Data set from J. Stephenson"
            self.shortname = "Steph."
            self.Stephenson()
            #if verbose: 
            print("Data set successfully loaded. {} trajectories.".format(self.size))
        else:
            print("There is an Error. You tried to load a data file of real distribution, but the file was not recognized.")
            
    def WD11(self):
        """ the data set is the Waszek and Deuss 2011 in the file WD11.dat """
        self.slices = ["PKIKP-PKiKP travel time residual", "turn lat",
                       "turn lon", "turn depth", "in lat", "in lon", "out lat", "out lon"]
        self.data = pd.read_table(self.filename, sep='\s+', names=self.slices, skiprows=10)#read_from_file(self.filename)
        if self.limited_nber_points[0]==True:
            print(self.limited_nber_points[1])
            self.data = self.data.iloc[:self.limited_nber_points[1]]
        self.size = self.data.shape[0]
        self.data_points = []
        for _, row in self.data.iterrows():
            ray = positions.Raypath()
            ray.add_b_t_point(positions.SeismoPoint(
                1. - row["turn depth"] / self.rICB, row["turn lat"], row["turn lon"]))
            in_point = positions.SeismoPoint(1., row["in lat"], row["in lon"])
            out_point = positions.SeismoPoint(
                1., row["out lat"], row["out lon"])
            ray.add_property({'in_point':in_point, 'out_point':out_point})
            ray.residual = row["PKIKP-PKiKP travel time residual"]
            self.data_points = np.append(self.data_points, ray)

    def Stephenson(self):  #the "turn depth" is actually the "turn radius" ! 
        self.slices = ["turn lat", "turn lon", "turn depth", "in lat", "in lon", 
                       "out lat", "out lon", "travel time residual relative to ak135"]
        nb_slices = [11,12,13,14,15,16,17,24]# [12,13,14,15,16,17,18,24]
        self.data = pd.read_table(self.filename, sep='\s+', header=None)[nb_slices]
        if self.limited_nber_points[0]==True:
            print(self.limited_nber_points[1])
            self.data = self.data.iloc[:self.limited_nber_points[1]]
        self.data.columns = self.slices
        self.size = self.data.shape[0]
        for _, row in self.data.iterrows():
            ray = positions.Raypath()
            ray.add_b_t_point(positions.SeismoPoint(
                row["turn depth"] / self.rICB, row["turn lat"], row["turn lon"]))
            in_point = positions.SeismoPoint(1., row["in lat"], row["in lon"])
            out_point = positions.SeismoPoint(
                1., row["out lat"], row["out lon"])
            ray.add_property({'in_point':in_point, 'out_point':out_point})
            ray.residual = row["travel time residual relative to ak135"]  # careful here with the names of the column for residual!
            self.data_points = np.append(self.data_points, ray)

    def real_residual(self):
        """ Extract the values of residuals from the data. """
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
        self.shortname = "equatorial"
        for x in np.linspace(-self.rICB, self.rICB, N):
            for y in np.linspace(-self.rICB, self.rICB, N):
                ray = positions.Raypath()
                ray.add_b_t_point(positions.CartesianPoint(x, y, 0.))
                if ray.bottom_turning_point.r <= self.rICB:
                    self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)

    def plot_c_vec(self, modelgeodyn, proxy=1, cm=plt.get_cmap('summer'), nameproxy=""):
        """ Plot contourf of the proxy + streamlines of the flow in meridional cross section.

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

    def plot_c(self, modelgeodyn, proxy=1, cm=plt.get_cmap('summer'), nameproxy=""):
        """ Plot contourf of the proxy in meridional cross section -- no stream lines.

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

        theta = np.linspace(0., 2 * np.pi, 1000)
        ax.plot(np.sin(theta), np.cos(theta), 'k', lw=3)
        ax.set_xlim([-1.1, 1.1])
        ax.set_ylim([-1.1, 1.1])

        cbar = plt.colorbar(sc)
        cbar.set_label(nameproxy)
        title = "Geodynamical model: {}".format(modelgeodyn.name)
        plt.title(title)
        plt.axis("off")


class RandomData(SeismicData):
    """ Random repartition of point, between depth 106 and 15km"""

    def __init__(self, N, rICB=1.):
        SeismicData.__init__(self)
        self.rICB = rICB
        self.N = N
        self.name = "Random repartition of data, between 0 and 100km depth"
        self.shortname = "random_0-100"
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
        self.shortname = "equatorialplane"
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


class Equator_upperpart(SeismicData):
    """ meshgrid for the uppermost part of IC (2-120km), at the equator."""

    def __init__(self, Nr, Np, rICB=1., d0=2., d1=120.):
        SeismicData.__init__(self)
        self.rICB = rICB
        self.N = Nr*Np
        self.Np = Np
        self.Nr = Nr
        self.name = "Meshgrid at the equator between 0 and 120km depth"
        self.shortname = "meshgrid"
        self.depth = [d0 / 1221., d1 / 1221.]
        self.theta = 0. # at the equator
        for depth in np.linspace(self.depth[0], self.depth[1], Nr):
            for phi in np.linspace(-180., 180., Np):
                ray = positions.Raypath()
                point = positions.SeismoPoint(self.rICB-depth, 0., phi)
                ray.add_b_t_point(point)
                self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)

    def mesh_RPProxy(self, proxy):
        r, t, p = self.extract_rtp("bottom_turning_point")
        #depth = (self.rICB - r)*self.rICB
        R = r.reshape(-1, self.Np)
        PHI = p.reshape(-1, self.Np)
        PROXY = proxy.reshape(-1, self.Np)
        return R, PHI, PROXY 


class PerfectSamplingSurface(SeismicData):

    def __init__(self, N, depth=0., rICB=1.):
        """ Grid of points partitioned at the surface (or at given depth under surface)
            
        :: arg depth:: depth in percentage below ICB.

        """
        SeismicData.__init__(self)
        self.rICB = rICB
        self.N = N
        self.name = "Perfect sampling at the surface"
        self.shortname = "surface"
        for t in np.linspace(-90, 90, N):
            for p in np.linspace(-180, 180, N):
                ray = positions.Raypath()
                ray.add_b_t_point(positions.SeismoPoint(rICB-rICB*depth, t, p))
                if ray.bottom_turning_point.r <= self.rICB:
                    self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)


    def mesh_TPProxy(self, proxy):
        r, t, p = self.extract_rtp("bottom_turning_point")
        #depth = (self.rICB - r)*self.rICB
        THETA = t.reshape(-1, self.N)
        PHI = p.reshape(-1, self.N)
        PROXY = proxy.reshape(-1, self.N)
        return THETA, PHI, PROXY 



# Joanne Stephenson data set: 
#12 - bottoming latitude
#13 - bottoming longitude
#14 - bottoming depth
#15 - IC entry lat 
#16 - IC entry lon
#17 - IC exit lat
#18 - IC exit lon
#25 - travel time residual relative to ak135



# This function is way to complicated, and reshape work well for what we need. Keep it here in case we actually have randomly distributed points on a grid. 
# def construct_meshgrid(x, y, z, info_x, info_y):
#     """ Construct a meshgrid based on x, y and z, to plot contourf.
#     x, y, z: 1D arrays, same size
#     x, y are the positions, and z the value that you want to plot.
#     info_x and info_y are the info to construct X and Y, using linspace.
#     return X, Y, Z that can be used directly for contour and contourf.
#     """
#     min_x, max_x, N_x = info_x
#     min_y, max_y, N_y = info_y
#     x1 = np.linspace(min_x, max_x, N_x)
#     y1 = np.linspace(min_y, max_y, N_y)
#     X, Y = np.meshgrid(x1, y1)
#     Z = -1. * np.ones_like(X)
#     print(x1.shape, y1.shape, X.shape, Y.shape, z.shape)
#     for it, z_element in enumerate(z):
#         ix = [i for i, j in enumerate(x1) if j == x[it]]
#         iy = [i for i, j in enumerate(y1) if j == y[it]]
#         Z[ix, iy] = z_element
#     return X, Y, Z
# 

