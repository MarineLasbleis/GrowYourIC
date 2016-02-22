#!/usr/local/bin/python
# Time-stamp: <2016-02-22 22:26:01 marine>
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


if __name__ == '__main__':


    #Initialize map
    m, fig = plot_data.setting_map()
    cm = plt.cm.get_cmap('RdYlBu')
    
    #seismic data set (from Lauren's file)
    data_points = read_write.read_from_file("results.dat", columns=[13, 14, 15], lines=-1)
    nlines, ncolumns = data_points.shape

    #translate it in the correct format
    dataset = []
    translation_dataset = []
    phi, theta, age = [], [], []
    print nlines, "points to write."
    for i in range(nlines):
        if i%100==0: print "Writing point ",  i, ". Coordinates: ", data_points[i, :]
        dataset.append(positions.Raypath("BT-point",
                                         1221.-data_points[i,2],
                                         data_points[i, 0],
                                         data_points[i, 1]))
        translation_dataset.append(geodynamic.Translation(positions.Point(1221.-data_points[i,2],
                                                          data_points[i, 0],
                                                          data_points[i, 1],
                                                          "seismo")))
        translation_dataset[i].analytical(2.01*1221./1., 1.)


        phi.append(translation_dataset[i].initial_position.phi)
        theta.append(translation_dataset[i].initial_position.theta)
        age.append(translation_dataset[i].exact_solution)

    
    x, y = m(phi, theta)
    m.scatter(x, y, c=age, zorder=10, cmap=plt.cm.RdYlGn)
    m.colorbar()
    plt.show()
