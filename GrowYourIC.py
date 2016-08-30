#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis

import numpy as np
import matplotlib.pyplot as plt #for figures

import positions
import geodyn, geodyn_trg, geodyn_static
import data

if __name__ == '__main__':
    ## Non-dimensional variables
    rICB = 1.
    age_ic = 1.
    omega = -0.5*np.pi # over write rotation rate. Rotation rates has to be in ]-np.pi, np.pi[
    # Define the velocity: 
    velocity_amplitude = 1.1
    velocity_center = [0., 100.]#center of the eastern hemisphere
    velocity = geodyn_trg.translation_velocity(velocity_center, velocity_amplitude)

    geodynModel = geodyn_trg.PureTranslation()
#    geodynModel = geodyn_trg.TranslationRotation()
    geodynModel = geodyn_trg.PureGrowth()
    geodynModel = geodyn_trg.TranslationGrowth()
#    geodynModel = geodyn_trg.TranslationGrowthRotation()
#    geodynModel = geodyn_static.Hemispheres()

    parameters = {'units': None, #non-dimensional 
                  'rICB': rICB, 
                  'tau_ic':age_ic,
                  'vt': velocity,
                  'exponent_growth': 0.5,
                  'omega': omega, 
                  'proxy_type': "age"}

    geodynModel.set_parameters(parameters)
    geodynModel.define_units()
    print geodynModel.__dict__

    ##  perfect sampling equator
    npoints = 20 #number of points in the x direction for the data set. 
    data_set = data.PerfectSamplingEquator(npoints)
    data_set.method = "bt_point"
    proxy = geodyn.evaluate_proxy(data_set, geodynModel)
    data_set.proxy = proxy #evaluate_proxy(data_set, geodynModel)
    data_set.plot_c_vec(geodynModel, proxy=proxy)
    #data_set.plot_scatter()
   # plt.show()

    data_set2 = data.PerfectSamplingEquatorRadial(100, 9)
    data_set2.method = "bt_point"
    proxy = geodyn.evaluate_proxy(data_set2, geodynModel)
    data_set2.proxy = proxy #evaluate_proxy(data_set, geodynModel)
    data_set2.radius_plot(geodynModel)

# #   
# # random data set
# data_set_random = data.RandomData(3000)
# data_set_random.method = "bt_point"
# proxy_random = geodyn.evaluate_proxy(data_set_random, geodynModel)
# data_set_random.proxy = proxy_random
# data_set_random.map_plot(geodynModel.name)
# data_set_random.phi_plot(geodynModel.name)
# data_set_random.distance_plot(geodynModel.name, positions.SeismoPoint(1., 0., -80.))
#  plt.show()
# 
    ## real data set
    data_set2 = data.SeismicFromFile("results.dat")
    data_set2.method = "bt_point"
    proxy2 = geodyn.evaluate_proxy(data_set2, geodynModel)
    data_set2.proxy = proxy2 #evaluate_proxy(data_set, geodynModel)
    data_set2.map_plot(geodynModel.name)
    data_set2.phi_plot(geodynModel.name)
    data_set2.distance_plot(geodynModel.name, positions.SeismoPoint(1., 0., -80.))
     #plt.show()
# 
# data_set2.proxy = data_set2.real_residual()
# data_set2.phi_plot(geodynModel.name)
# data_set2.map_plot(geodynModel.name)
# data_set2.distance_plot(geodynModel.name, positions.SeismoPoint(1., 0., -80.))
# 
#    ## real data set, average over raypath
# data_set3 = data.SeismicFromFile("results.dat")
# data_set3.method = "raypath"
# geodynModel.evaluation= "0"
# data_set3.NpointsRaypath = 20 
# proxy3 = geodyn.evaluate_proxy(data_set3, geodynModel)
# data_set3.proxy = proxy3 #evaluate_proxy(data_set, geodynModel)
# data_set3.map_plot(geodynModel.name)
# data_set3.phi_plot(geodynModel.name)
# data_set3.distance_plot(geodynModel.name, positions.SeismoPoint(1., 0., -80.))
# 
    plt.show()
