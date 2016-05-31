#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis

import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math

import positions
import geodynamic
import plot_data
import data

## Choose a model in the list: 
    # geodynamic.PureTranslation(velocity_translation)
    # geodynamic.TranslationRotation(velocity_translation, omega)



## Choose a data set:
    # data.SeismicFromFile(filename)
        # plot function is : 
    # data.PerfectSamplingEquator(numbers_of_points)
        # plot function is : 

## Choose a visualization tool:
    # surface plot
    # meridional section
    # equatorial section
    
    
## TODO : adim ray in positions.py


if __name__ == '__main__':
    
    # rICB = 1221.   #km
    # age_ic = 1.e9  # years
    # velocity = rICB/90e6*2. *np.array([1., 0., 0.])# translation velocity. km/years
    # omega = 0.00001* 2.*np.pi #rad/years

    # # Non-dimensionalisation of the variables
    # velocity = velocity*age_ic/rICB
    # omega = omega*age_ic
    rICB = 1.
    age_ic = 1.

    velocity = [1, 0., 0.]
    omega = -0.5*np.pi # over write rotation rate. Rotation rates has to be in ]-np.pi, np.pi[

    print "velocity: ", velocity
    print "omega: ", omega, "en radians"
    print "age", age_ic, ", rICB: ", rICB



    geodynModel = geodynamic.PureTranslation()
    geodynModel = geodynamic.TranslationRotation()
    geodynModel = geodynamic.PureGrowth()
    geodynModel = geodynamic.TranslationGrowth()
    #geodynModel = geodynamic.TranslationGrowthRotation()
    geodynModel.set_vt(velocity)
    geodynModel.set_tauIC(age_ic) # made dimensionless by using age_ic
    geodynModel.set_exponent_growth(0.3)
    geodynModel.set_rICB(rICB)
    geodynModel.set_rotation(omega)
    
    # geodynModel.plot_equatorial(-1, 1)

    ##  perfect sampling equator
    npoints = 50 #number of points in the x direction for the data set. 
    data_set = data.PerfectSamplingEquator(npoints, rICB = 1.)
    #this data set is already dimensionless!
    data_set.method = "bt_point"
    proxy = geodynamic.evaluate_proxy(data_set, geodynModel)
    data_set.proxy = proxy #evaluate_proxy(data_set, geodynModel)
    data_set.plot_c_vec(geodynModel)
    data_set.plot_scatter()
#
    plt.show()
#
#    ## real data set
#    data_set2 = data.SeismicFromFile("results.dat")
#    data_set2.method = "bt_point"
#    proxy2 = geodynamic.evaluate_proxy(data_set2, geodynModel)
#    data_set2.proxy = proxy2 #evaluate_proxy(data_set, geodynModel)
#    data_set2.map_plot()
#    data_set2.phi_plot()
#
#    ## real data set
#    data_set3 = data.SeismicFromFile("results.dat")
#    data_set3.method = "raypath"
#    data_set3.NpointsRaypath = 20 
#    proxy3 = geodynamic.evaluate_proxy(data_set3, geodynModel)
#    data_set3.proxy = proxy3 #evaluate_proxy(data_set, geodynModel)
#    data_set3.map_plot()
#    data_set3.phi_plot()
#
#    plt.show()
