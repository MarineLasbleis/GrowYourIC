#!/usr/local/bin/python
# Time-stamp: <2016-03-08 07:42:35 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

#personal routines
import positions
import plot_data
import read_write
import geodynamic



def Translation_analytic_BT(velocity_translation = 2.01*1221./1.e6, direction=positions.CartesianPoint(1,0,0)):

    #Initialize map
    m, fig = plot_data.setting_map()
    cm = plt.cm.get_cmap('RdYlBu')
    
    #seismic data set (from Lauren's file)
    data_points = read_write.read_from_file("results.dat", slices=["PKIKP-PKiKP travel time residual", "turn lat", "turn lon", "turn depth"])
    nlines, ncolumns = data_points.shape
    print data_points.info
    #translate it in the correct format
    dataset = []
    translation_dataset = []
    phi, theta, age, dt = [], [], [], []
    print nlines, "points to write."
    for i in range(nlines):
        if i%100==0: #print every 100 values
            print "Writing point ",  i, ". Coordinates: ", data_points.ix[i]

        # because it's analytical solution, we don't need to define a grid. model is calulated exactly. 
        translation_dataset.append(geodynamic.Translation(positions.SeismoPoint(1221.-data_points.ix[i,"turn depth"],
                                                                                 data_points.ix[i, "turn lat"],
                                                                                 data_points.ix[i, "turn lon"])))
        translation_dataset[i].analytical(velocity_translation, direction)


        phi.append(translation_dataset[i].initial_position.phi)
        theta.append(translation_dataset[i].initial_position.theta)
        age.append(translation_dataset[i].exact_solution)

    x, y = m(phi, theta)
    m.scatter(x, y, c=age, zorder=10, cmap=plt.cm.RdYlGn)
    m.colorbar()

    fig, ax = plt.subplots()
    ax.plot(phi, age/max(age), ".")
    dt = data_points["PKIKP-PKiKP travel time residual"]
    print dt.shape
    ax.plot(phi, dt,".r")

    plt.show()



def Translation_analytical_raypath(velocity_translation = 2.01*1221./1.e6, direction=positions.CartesianPoint(1,0,0)):

    #Initialize map
    m, fig = plot_data.setting_map()
    cm = plt.cm.get_cmap('RdYlBu')

    #seismic data set (from Lauren's file)
    data_points = read_write.read_from_file("results.dat", slices=["PKIKP-PKiKP travel time residual", "turn lat", "turn lon", "turn depth", "in lat", "in lon", "out lat", "out lon",])
    nlines, ncolumns = data_points.shape
    print data_points.info

    #translate it in the correct format
    dataset = []
    translation_dataset = []
    phi, theta, age, dt = [], [], [], []
    print nlines, "points to write."
    for i in range(nlines):
        if i%100==0: #print every 100 values
            print "Writing point ",  i, ". Coordinates: ", data_points.ix[i]
        # because it's analytical solution, we don't need to define a grid. model is calulated exactly.
        translation_dataset.append(geodynamic.Translation(positions.SeismoPoint(1221.-data_points.ix[i,"turn depth"],
                                                                                 data_points.ix[i, "turn lat"],
                                                                                 data_points.ix[i, "turn lon"])))
        translation_dataset[i].analytical(velocity_translation, direction)

        phi.append(translation_dataset[i].initial_position.phi)
        theta.append(translation_dataset[i].initial_position.theta)
        age.append(translation_dataset[i].exact_solution)

    x, y = m(phi, theta)
    m.scatter(x, y, c=age, zorder=10, cmap=plt.cm.RdYlGn)
    m.colorbar()

    fig, ax = plt.subplots()
    ax.plot(phi, age/max(age), ".")
    dt = data_points["PKIKP-PKiKP travel time residual"]
    print dt.shape
    ax.plot(phi, dt,".r")

    plt.show()

                                

if __name__ == '__main__':

    age_IC = 1.e6 #in years
    velocity_translation = 2.01*1221./1.e6 #in km/years
    direction = positions.SeismoPoint(1, 0, 80)

    # Translation_analytic_BT(velocity_translation, direction)

    Translation_analytical_raypath()
