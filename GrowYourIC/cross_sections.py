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


if __name__ == '__main__':

    
    rICB = 1.
    age_ic = 1.
    velocity = [.5, 0., 0.]
    omega = -0.5*np.pi # over write rotation rate. Rotation rates has to be in ]-np.pi, np.pi[


    models = [geodyn_trg.PureTranslation(), geodyn_trg.TranslationRotation(), geodyn_trg.PureGrowth(), geodyn_trg.TranslationGrowth(), geodyn_trg.TranslationGrowthRotation()] 

    models = [geodyn_trg.TranslationGrowth(), geodyn_trg.TranslationGrowthRotation()]
    geodynModel = geodyn_trg.PureTranslation()
    geodynModel = geodyn_trg.TranslationRotation()
    geodynModel = geodyn_trg.PureGrowth()
    geodynModel = geodyn_trg.TranslationGrowth()
#    geodynModel = geodyn_trg.TranslationGrowthRotation()
#    geodynModel = geodyn_static.Hemispheres()

    for geodynModel in models:
        
        parameters = { 'rICB': rICB, 
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
    plt.show()

