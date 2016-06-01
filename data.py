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
    
    def extract_xyz(self, type_of_point):
        assert self.size, 'data_points is probably empty' # TO DO : raise exceptions ins    tead of using assert
        x, y, z = np.empty([self.size, 1]), np.empty([self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            point = getattr(ray, type_of_point)
            x[i] = point.x
            y[i] = point.y
            z[i] = point.z
        return x, y, z 

    def extract_rtp(self, type_of_point):
        """Extract the radius, theta (latitute), phi (longitude) for a serie of points 
    
        """
        assert self.size, 'data_points is probably empty' # TO DO : raise exceptions ins    tead of using asse    rt
        r, theta, phi = np.empty([self.size, 1]), np.empty([self.size, 1]), np.empty([self.size, 1])
        for i, ray in enumerate(self.data_points):
            point = getattr(ray, type_of_point)
            r[i] = point.r
            theta[i] = point.theta
            phi[i] = point.phi
        return r, theta, phi

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

    def map_plot(self, geodyn_model=''):
        """ plot data on a map."""
        # need to check which data exist, if raypath, BT point, in-out point, etc. 
        # and it should plot as much data as possible (only BT if only exist, but raypath also)
        # should also ask for which data set you want to plot (or plot all of them? dt from data, and results if they exist?)

        m, fig = plot_data.setting_map()
        cm = plt.cm.get_cmap('RdYlBu')
        
        r, theta, phi = self.extract_btpoints()
        r, theta, phi = self.extract_rtp("bottom_turning_point")
        x, y = m(phi, theta)
        proxy = np.array([self.proxy]).T.astype(float)
        sc = m.scatter(x, y, c=proxy, zorder=10, cmap=cm)
        
        # TO DO : make a function to plot great circles correctly!
        #r1, theta1, phi1 = self.extract_in()
        #r2, theta2, phi2 = self.extract_out()
        #for i, t in enumerate(theta1):
        #    z, w = m.gcpoints(phi1[i], theta1[i], phi2[i], theta2[i], 200)#
        #    m.plot(z, w, zorder=5, c="black")
        #    m.drawgreatcircle(phi1[i], theta1[i], phi2[i], theta2[i], zorder=5, c="black")
        title = "Dataset: {},\n geodynamic model: {}".format(self.name, geodyn_model)
        plt.title(title)
        plt.colorbar(sc)
        #plt.show()


    def phi_plot(self, geodyn_model=''):
        """ Plot proxy as function of longitude """

        fig, ax = plt.subplots()
        r, theta, phi = self.extract_rtp("bottom_turning_point")
        ax.plot(phi, self.proxy, '.')
        title = "Dataset: {},\n geodynamic model: {}".format(self.name, geodyn_model)
        plt.title(title)
        plt.xlabel("longitude of bottom turning point")
        plt.ylabel("proxy")
        #plt.show()

    def adimension(self, scale):
        """ coordinates of points in the data set are made dimensionless using the value of scale """

        if self.size == None :
            raise IndexError("The data set is empty!")
        elif self.size <= 0:
            raise IndexError("The data set is empty!")
        else:
            for i, ray in self.data_points:
                self.data_points[i].adim(scale)


class SeismicFromFile(SeismicData):

    def __init__(self, filename="results.dat", RICB=1221.):
        
        SeismicData.__init__(self)

        self.name = "Data set from Lauren's file"
        #seismic data set (from Lauren's file)
        self.filename = filename
        self.slices = ["PKIKP-PKiKP travel time residual", "turn lat", "turn lon", "turn depth", "in lat", "in lon", "out lat", "out lon"]
        self.data = read_from_file(filename, slices=self.slices)
        self.size = self.data.shape[0]

        self.data_points = []

        for i, row in self.data.iterrows():
            ray = positions.Raypath()
            ray.add_b_t_point(positions.SeismoPoint(1.-row["turn depth"]/RICB, row["turn lat"], row["turn lon"]))
            in_Point = positions.SeismoPoint(1., row["in lat"], row["in lon"])
            out_Point = positions.SeismoPoint(1., row["out lat"], row["out lon"])
            ray.add_in_out(in_Point, out_Point)
            self.data_points = np.append(self.data_points, ray)
            #self.data_points.append(ray)
        assert(self.size == len(self.data_points))

class PerfectSamplingEquator(SeismicData):
    
    def __init__(self, N, rICB = 1221.):
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

    def plot_scatter(self):
        fig, ax = plt.subplots()
        x = np.linspace(-self.rICB, self.rICB , 100)
        ax.plot(x, np.sqrt(self.rICB**2-x**2), 'k')
        ax.plot(x, -np.sqrt(self.rICB**2-x**2), 'k')
        if hasattr(self, "proxy"):
            proxy = np.array([self.proxy]).T.astype(float)
        else :
            proxy = 1.
        x, y, z = self.extract_xyz("bottom_turning_point")
        sc = ax.scatter(x, y, c=proxy)
        plt.axis("equal")
        plt.colorbar(sc)
        plt.show()

    def plot_contourf(self):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        x = np.linspace(-self.rICB, self.rICB , 100)
        ax.plot(x, np.sqrt(self.rICB**2-x**2), 'k')
        ax.plot(x, -np.sqrt(self.rICB**2-x**2), 'k')
        if hasattr(self, "proxy"):
            proxy = self.proxy
        else :
            proxy = 1.
        x1 = np.linspace(-self.rICB, self.rICB , self.N)
        y1 = np.linspace(-self.rICB, self.rICB , self.N)
        X, Y = np.meshgrid(x1, y1)
        Z = -1.*np.ones_like(X)
        x, y, z = self.extract_xyz("bottom_turning_point")
        for it, pro in enumerate(proxy):
            ix = [i for i, j in enumerate(x1) if j == x[it]]
            iy = [i for i, j in enumerate(y1) if j == y[it]]
            Z[ix, iy] = pro
        Z = np.ma.array(Z, mask = Z==-1)    
        sc = ax.contourf(Y, X, Z, 100)
        
        plt.colorbar(sc)
        #plt.show()


    def plot_c_vec(self, modelgeodyn):

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        x = np.linspace(-self.rICB, self.rICB , 100)
        ax.plot(x, np.sqrt(self.rICB**2-x**2), 'k')
        ax.plot(x, -np.sqrt(self.rICB**2-x**2), 'k')
        if hasattr(self, "proxy"):
            proxy = self.proxy
        else :
            proxy = 1.
        x1 = np.linspace(-self.rICB, self.rICB , self.N)
        y1 = np.linspace(-self.rICB, self.rICB , self.N)
        X, Y = np.meshgrid(x1, y1)
        Z = -1.*np.ones_like(X)
        x, y, z = self.extract_xyz("bottom_turning_point")
        for it, pro in enumerate(proxy):
            ix = [i for i, j in enumerate(x1) if j == x[it]]
            iy = [i for i, j in enumerate(y1) if j == y[it]]
            Z[ix, iy] = pro
        mask_Z = Z==-1    
        Z = np.ma.array(Z, mask = mask_Z)
        sc = ax.contourf(Y, X, Z, 100)
        
        Vx, Vy = np.empty((self.N, self.N)), np.empty((self.N, self.N))
        for ix, xi in enumerate(x1):
            for iy, yi in enumerate(y1):
                velocity = modelgeodyn.velocity(modelgeodyn.tau_ic, [X[ix, iy], Y[ix, iy], 0.])
                Vx[ix, iy] = velocity[0]
                Vy[ix, iy] = velocity[1]
        Vx = np.ma.array(Vx, mask=mask_Z)
        Vy = np.ma.array(Vy, mask=mask_Z)
        #ax.quiver(X, Y, Vx, Vy)
        ax.streamplot(X,Y,Vx,Vy, color='black', arrowstyle = '->')
        plt.colorbar(sc)
        title = "Geodynamical model: {}".format(modelgeodyn.name)
        plt.title(title)
        plt.axis("off")
        #plt.show()

class RandomData(SeismicData):

    def __init__(self, N, rICB = 1.):
        SeismicData.__init__(self)
        self.rICB = rICB
        self.N = N
        self.name = "Random repartition of data, between 0 and 100km depth"
        self.random_method = "uniform"
        self.depth = [0., 100./1221.]
        
        for i in range(N):
            ray = positions.Raypath()
            ray.add_b_t_point(positions.RandomPoint(self.random_method, self.depth, rICB))
            self.data_points = np.append(self.data_points, ray)
        self.size = len(self.data_points)    



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
