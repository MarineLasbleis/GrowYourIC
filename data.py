#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import pandas as pd

# personal routines
import positions
import geodynamic
import plot_data

RICB = 1221.

def read_from_file(filename, names=["station", "PKIKP-PKiKP travel time residual", "zeta", "epicentral distance", "station lat", "stat    ion lon", "event lat", "event lon", "event depth", "in lat", "in lon", "out lat", "out lon", "turn lat", "turn lon", "turn depth", "in    ner core travel time", "PKIKP/PKiKP amplitude ratio"], slices="all"):
    """ read seismic data repartition
    
    input parameters:
    - filename: name of the data file
    - names: names of the columns for the data set
    - slices: names of columns for the output.
    output:
    - data : pandas DataFrame with all the datas. Columns name are indicated by the variable "names".
    """
    df = pd.read_table(filename, sep=' ', names=names, skiprows=0)
    if slices != "all":
        df = df[slices]
    return df


class SeismicData():
    """ Class for seismic data """
    
    def __init__(self):
        self.data_points = [] 
        self.size = None
    
    def __getitem__(self, key):
        return self.data_points[key]
    
    def extract_btpoints(self):
        assert self.size, 'data_points is probably empty' # TO DO : raise exceptions instead of using assert
        # need to also assert that bottom_turning_point exist or can be calculated!
        r, theta, phi = np.empty([self.size, 1]), np.empty([self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            r[i] = ray.bottom_turning_point.r
            theta[i] = ray.bottom_turning_point.theta
            phi[i] = ray.bottom_turning_point.phi
        return r, theta, phi 

    def extract_in(self):
        assert self.size, 'data_points is probably empty' # TO DO : raise exceptions instead of using assert
        # need to also assert that bottom_turning_point exist or can be calculated!
        r, theta, phi = np.empty([self.size, 1]), np.empty([self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            r[i] = ray.in_point.r
            theta[i] = ray.in_point.theta
            phi[i] = ray.in_point.phi
        return r, theta, phi

    def extract_out(self):
        assert self.size, 'data_points is probably empty' # TO DO : raise exceptions instead of using assert
        # need to also assert that bottom_turning_point exist or can be calculated!
        r, theta, phi = np.empty([self.size, 1]), np.empty([self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            r[i] = ray.out_point.r
            theta[i] = ray.out_point.theta
            phi[i] = ray.out_point.phi
        return r, theta, phi

    def translation_BT(self, velocity, direction):
        assert self.size, 'data_points is probably empty' # TO DO : raise exceptions instead of using assert
        # need to also assert that bottom_turning_point exist or can be calculated!
        self.translation = np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            self.translation[i] = geodynamic.exact_translation(ray.bottom_turning_point, velocity, direction)
    
    def translation_raypath(self, velocity, direction, N=10):
        """ 
        N : number of points in the trajectory
        """
        # need to check in raypath exist. In case it does not, need to apply one of the method for raypath
        self.translation = np.zeros([self.size, 1])
        for i, ray in enumerate(self.data_points):
            #assuming raypath does not exist, so need to be calculated (but would be faster if it checks before and do not run this if it's not needed)
            self.data_points[i].straigth_in_out(N)
            raypath = ray.points #raypath is a np array, each elements being one point.
            total_translation = 0. 
            for j, points in enumerate(raypath):
                _translation = geodynamic.exact_translation(points, velocity, direction)
                total_translation += _translation

            #total_translation /= N #TO DO : add the weigth by the distance (each points may be at different distances from each others?)
            self.translation[i] = total_translation /float(N)

    def map_plot(self):
        """ plot data on a map."""
        #need to check which data exist, if raypath, BT point, in-out point, etc. 
        # and it should plot as much data as possible (only BT if only exist, but raypath also)
        # should also ask for which data set you want to plot (or plot all of them? dt from data, and results if they exist?)

        m, fig = plot_data.setting_map()
        cm = plt.cm.get_cmap('RdYlBu')
        
        r, theta, phi = self.extract_btpoints()

        x, y = m(phi, theta)
        m.scatter(x, y, c=self.translation, zorder=10, cmap=cm)
        
        # TO DO : make a function to plot great circles correctly!
        #r1, theta1, phi1 = self.extract_in()
        #r2, theta2, phi2 = self.extract_out()
        #for i, t in enumerate(theta1):
        #    z, w = m.gcpoints(phi1[i], theta1[i], phi2[i], theta2[i], 200)#
        #    m.plot(z, w, zorder=5, c="black")
        #    m.drawgreatcircle(phi1[i], theta1[i], phi2[i], theta2[i], zorder=5, c="black")

        #plt.show()




class SeismicFromFile(SeismicData):

    def __init__(self, filename="results.dat"):
        
        SeismicData.__init__(self)

        #seismic data set (from Lauren's file)
        self.filename = filename
        self.slices = ["PKIKP-PKiKP travel time residual", "turn lat", "turn lon", "turn depth", "in lat", "in lon", "out lat", "out lon"]
        self.data = read_from_file(filename, slices=self.slices)
        self.size = self.data.shape[0]

        self.data_points = []

        for i, row in self.data.iterrows():
            ray = positions.Raypath()
            ray.add_b_t_point(positions.SeismoPoint(RICB-row["turn depth"], row["turn lat"], row["turn lon"]))
            in_Point = positions.SeismoPoint(RICB, row["in lat"], row["in lon"])
            out_Point = positions.SeismoPoint(RICB, row["out lat"], row["out lon"])
            ray.add_in_out(in_Point, out_Point)
            self.data_points.append(ray)



class RandomData(SeismicData):
    pass


if __name__ == '__main__':

    
    f = SeismicFromFile()
    #print f.data_points
    f.translation_BT(1e6, positions.SeismoPoint(1,0,0))
    f.map_plot()

    g = SeismicFromFile()
    g.translation_raypath(1e6, positions.SeismoPoint(1,0,0))
    g.map_plot()

    fig, ax = plt.subplots()
    ax.plot(g.translation, f.translation, '.')

    plt.show()
