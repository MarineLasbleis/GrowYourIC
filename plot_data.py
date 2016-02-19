#!/usr/local/bin/python
# Time-stamp: <2016-02-19 16:41:57 marine>
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

    data = read_write.read_from_file('results.dat', columns=[10, 11, 12, 13], lines=[1,2,3])

    fig =plt.figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='moll',lon_0=0.,resolution='c')
    m.drawcoastlines(linewidth=0.25)
    #m.drawcountries(linewidth=0.25)
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))
    m.drawmapboundary(fill_color='#99ffff')


    
    
    plt.show()
