#!usr/local/bin/python
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
import data


def test_Lauren():
    f = data.SeismicFromFile()
    f.translation_BT(1e6, positions.SeismoPoint(1,0,0))
    f.map_plot()
    f.translation_raypath(1e6, positions.SeismoPoint(1,0,0))
    f.map_plot()
     
    plt.show()


if __name__ == '__main__':
    
    f = data.PerfectSamplingEquator(20)
    f.plot()
