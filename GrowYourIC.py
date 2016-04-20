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


## Choose a data set:
    # data.SeismicFromFile(filename)
        # plot function is : 
    # data.PerfectSamplingEquator(numbers_of_points)
        # plot function is : 


def evaluate_proxy(dataset, method):
    """ evaluate the value of the proxy on all the points of the data set, using the choosen geodynamical method 
    
    dataset : a data.SeismicData object
    method : a geodynamic.ModelGeodynamic object
    """
    # TO DO : choose evaluate on raypath or on BT point
    time = np.empty_like(dataset.data_points)
    for i, ray in enumerate(dataset.data_points):
        if dataset.method == "bt_point":
            point = ray.bottom_turning_point
            time[i] = evaluate_singlepoint(point, method)[0]
        elif dataset.method == "raypath":
            N = dataset.NpointsRaypath
            dataset.data_points[i].straigth_in_out(N)
            raypath = ray.points
            total_proxy = 0.
            for j, point in enumerate(raypath):
                _proxy = evaluate_singlepoint(point, method)[0]
                total_proxy += _proxy
            time[i] = total_proxy / float(N)    

    return time

def evaluate_singlepoint(point, method):
    """ evaluate the proxy on a single data.Point instance, using the choosen method."""
    x, y, z = point.x, point.y, point.z
    time = method.find_time_beforex0([x, y, z], method.tau_ic, method.tau_ic)
    return method.tau_ic-time







if __name__ == '__main__':
    
    rICB = 1221.
    age_ic = 1.e9 
    velocity = 2. *rICB / 180e3 *np.array([0., 0., 1.])# translation velocity
    

    geodynModel = geodynamic.PureTranslation(velocity)
    geodynModel.set_tauIC(age_ic)


    ##  perfect sampling equator
    npoints = 10 #number of points in the x direction for the data set. 
    data_set = data.PerfectSamplingEquator(npoints)
    data_set.method = "bt_point"
    proxy = evaluate_proxy(data_set, geodynModel)
    data_set.proxy = proxy #evaluate_proxy(data_set, geodynModel)
    data_set.plot_contourf()

    ## real data set
    data_set2 = data.SeismicFromFile("results.dat")
    data_set2.method = "bt_point"
    proxy2 = evaluate_proxy(data_set2, geodynModel)
    data_set2.proxy = proxy2 #evaluate_proxy(data_set, geodynModel)
    data_set2.map_plot()

    ## real data set
    data_set3 = data.SeismicFromFile("results.dat")
    data_set3.method = "raypath"
    data_set3.NpointsRaypath = 3 
    proxy3 = evaluate_proxy(data_set3, geodynModel)
    data_set3.proxy = proxy3 #evaluate_proxy(data_set, geodynModel)
    data_set3.map_plot()
