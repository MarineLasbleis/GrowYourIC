#!/usr/local/bin/python
# Time-stamp: <2016-02-19 22:52:54 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis


import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

#personal routines
import positions
import read_write


# def 


if __name__ == '__main__':


    # columns=[9, 10, 11, 12] in, out
    # columns=[14, 15, 16] bottom turning point (16 is depth, no radius!)
    data = read_write.read_from_file('results.dat', columns=[13, 14, 15], lines=-1)
    nlines, ncolumns = data.shape
    fig =plt.figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='moll',lon_0=0.,resolution='c')
    m.drawcoastlines(linewidth=0.25)
    #m.drawcountries(linewidth=0.25)
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))
    m.drawmapboundary(fill_color='#99ffff')


    dataset = []
    print nlines
    for i in range(nlines):
        if data[i,1]<0:
            print i, data[i, :]
        dataset.append(positions.Raypath("BT-point", 1221-data[i,2], data[i, 0],  data[i, 1]))
        x, y =m(dataset[i].bottom_turning_point.phi, dataset[i].bottom_turning_point.theta)
        m.scatter(x,y, zorder=10)
        #print i, dataset[i].bottom_turning_point.theta, dataset[i].bottom_turning_point.phi

    
    
    
    plt.show()