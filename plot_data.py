#!/usr/local/bin/python
# Time-stamp: <2016-02-22 21:53:55 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis


import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

#personal routines
import positions


def setting_map():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m = Basemap(projection='moll',lon_0=0.,resolution='c')
    m.drawcoastlines(linewidth=0.25)
    #m.drawcountries(linewidth=0.25)
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))
    m.drawmapboundary(fill_color='#99ffff')
    return m, fig

