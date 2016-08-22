#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis

import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math

import positions
import geodyn, geodyn_trg, geodyn_static
import plot_data
import data

## Choose a model in the list: 
    #geodyn_trg.PureTranslation()
    #geodyn_trg.TranslationRotation()
    #geodyn_trg.PureGrowth()
    #geodyn_trg.TranslationGrowth()
    #geodyn_trg.TranslationGrowthRotation()
    #geodyn_static.Hemispheres()

## Choose a proxy type:
    # age
    # position (! will give a Point instance)
    # phi
    # theta
    # growth rate

## set the parameters for the model : geodynModel.set_parameters(parameters)

## Choose a data set:
    # data.SeismicFromFile(filename)
        # plot function is : 
    # data.PerfectSamplingEquator(numbers_of_points)
        # plot function is : 
    # data.RandomData(numbers_of_points)
        # plot function is :

## Choose a visualization tool:
    # surface plot
    # meridional section
    # equatorial section
    
    
## TODO : change the figure to choose the values : that's difficult, because the function took the function f which is already the difference between the two functions (and not both functions)
## TODO : change the way it chooses the interval (remember the last choice and try it) : I added a way to try different intervals. But need improvements.
## TODO : see if we can add models of proxy (artificial mapping of proxy in the IC)
## TODO : proxy = growth rate at the cristallization time?
## TODO : proxy = position at the ICB at cristallisation time? (longitude?)
## TODO : variable growth rate (square root of time)
## TODO : effective_growth_rate : verify the function! 

if __name__ == '__main__':
    
    # rICB = 1221.   #km
    # age_ic = 1.e9  # years
    # velocity = rICB/90e6*2. *np.array([1., 0., 0.])# translation velocity. km/years
    # omega = 0.00001* 2.*np.pi #rad/years

    # # Non-dimensionalisation of the variables
    # velocity = velocity*age_ic/rICB
    # omega = omega*age_ic


    ## Non-dimensional variables

    rICB = 1.
    age_ic = 1.
    velocity = [3., 0., 0.]
    omega = -0.5*np.pi # over write rotation rate. Rotation rates has to be in ]-np.pi, np.pi[

    print "velocity: ", velocity
    print "omega: ", omega, "en radians"
    print "age", age_ic, ", rICB: ", rICB

# Define the velocity: 
    velocity_amplitude = 3.
    velocity_center = [0., 100.]#center of the eastern hemisphere
    velocity = geodyn_trg.translation_velocity(velocity_center, velocity_amplitude)

    geodynModel = geodyn_trg.PureTranslation()
#    geodynModel = geodyn_trg.TranslationRotation()
#    geodynModel = geodyn_trg.PureGrowth()
    geodynModel = geodyn_trg.TranslationGrowth()
#    geodynModel = geodyn_trg.TranslationGrowthRotation()
#    geodynModel = geodyn_static.Hemispheres()

    parameters = {'units':None, #non-dimensional 
                  'rICB': rICB, 
                  'tau_ic':age_ic,
                  'vt': velocity,
                  'exponent_growth': 0.3,
                  'omega': omega, 
                  'proxy_type': "age"}

    geodynModel.set_parameters(parameters)

    ##  perfect sampling equator
    npoints = 20 #number of points in the x direction for the data set. 
    data_set = data.PerfectSamplingEquator(npoints, rICB = 1.)
    data_set.method = "bt_point"
    proxy = geodyn.evaluate_proxy(data_set, geodynModel)
    data_set.proxy = proxy #evaluate_proxy(data_set, geodynModel)
    data_set.plot_c_vec(geodynModel)
    #data_set.plot_scatter()
   # plt.show()

   
    # random data set
    data_set_random = data.RandomData(300)
    data_set_random.method = "bt_point"
    proxy_random = geodyn.evaluate_proxy(data_set_random, geodynModel)
    data_set_random.proxy = proxy_random
    data_set_random.map_plot(geodynModel.name)
    data_set_random.phi_plot(geodynModel.name)
    #plt.show()
#
    ## real data set
    data_set2 = data.SeismicFromFile("results.dat")
    data_set2.method = "bt_point"
    proxy2 = geodyn.evaluate_proxy(data_set2, geodynModel)
    data_set2.proxy = proxy2 #evaluate_proxy(data_set, geodynModel)
    data_set2.map_plot(geodynModel.name)
    data_set2.phi_plot(geodynModel.name)
     #plt.show()
# # 
# #  #     ## real data set, average over raypath
# #     data_set3 = data.SeismicFromFile("results.dat")
# #     data_set3.method = "raypath"
# #     data_set3.NpointsRaypath = 20 
# #     proxy3 = geodyn.evaluate_proxy(data_set3, geodynModel)
# #     data_set3.proxy = proxy3 #evaluate_proxy(data_set, geodynModel)
# #     data_set3.map_plot(geodynModel.name)
# #     data_set3.phi_plot(geodynModel.name)
# # # #
    plt.show()
