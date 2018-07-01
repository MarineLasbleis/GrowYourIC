import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math

from GrowYourIC import positions, geodyn, geodyn_trg, geodyn_static, plot_data, data

#plt.rcParams['figure.figsize'] = (8.0, 3.0) #size of figures
cm = plt.cm.get_cmap('viridis')

#data_set = data.SeismicFromFile("./GrowYourIC/data/WD11.dat", N=20)
data_set = data.SeismicFromFile("./GrowYourIC/data/DF_sample_ksi_sorted.dat", N=500)


r, t, p = data_set.extract_rtp("bottom_turning_point")

## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c="black",s=20, zorder=10, cmap=cm, edgecolors='none')


geodynModel = geodyn_trg.TranslationGrowthRotation()

parameters = {'units': None,
              'rICB': 1., 
              'tau_ic':1.,
              'vt': [2.5,0,0],  #has to be 3-D velocity
              'exponent_growth': 0.,
              'omega': 0.,
              'proxy_type': "age"}
geodynModel.set_parameters(parameters)
geodynModel.define_units()


npoints = 50 #number of points in the x direction for the data set. 
data_set_2 = data.PerfectSamplingEquator(npoints, rICB = 1.)
data_set_2.method = "bt_point"
proxy = geodyn.evaluate_proxy(data_set_2, geodynModel, proxy_type="age", verbose = False)
data_set_2.plot_c_vec(geodynModel, proxy=proxy, cm=cm, nameproxy="age (Myears)")



proxy_name = "age (Myears)"


data_set.method = "raypath"
data_set.NpointsRaypath = 10
proxy_real = geodyn.evaluate_proxy(data_set, geodynModel, proxy_type="age", verbose=False)

## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy_real,s=8, zorder=10, cmap=cm, edgecolors='none')
plt.title("Dataset: {},\n geodynamic model: {}".format(data_set.name, geodynModel.name))
cbar = plt.colorbar(sc)
cbar.set_label(proxy_name)





npoints = 60 #number of points in the x direction for the data set. 
data_set = data.PerfectSamplingEquator(npoints, rICB = 1.)
data_set.method = "bt_point"


Radial = geodyn_static.Radial_sym()
proxy = geodyn.evaluate_proxy(data_set, Radial, proxy_type="radius", verbose = False)
data_set.plot_c_vec(Radial, proxy=proxy, cm=cm, nameproxy="age (Myears)")

Radial = geodyn_static.Innermost_IC(0.5)
proxy = geodyn.evaluate_proxy(data_set, Radial, proxy_type="radius", verbose = False)
data_set.plot_c_vec(Radial, proxy=proxy, cm=cm, nameproxy="age (Myears)")

plt.show()

data_set = data.SeismicFromFile("./GrowYourIC/data/DF_sample_ksi_sorted.dat")
r, t, p = data_set.extract_rtp("bottom_turning_point")


data_set.method = "bt_point"
proxy = geodyn.evaluate_proxy(data_set, Radial, proxy_type="radius", verbose = False)
## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy+1,s=8, zorder=10, cmap=cm, edgecolors='none')
plt.title("Dataset: {},\n geodynamic model: {}".format(data_set.name, Radial.name))
#cbar = plt.colorbar(sc)


plt.show()

data_set.method = "raypath"
data_set.NpointsRaypath = 20
proxy = geodyn.evaluate_proxy(data_set, Radial, proxy_type="radius", verbose = False)
## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy,s=8, zorder=10, cmap=cm, edgecolors='none')
plt.title("Dataset: {},\n geodynamic model: {}".format(data_set.name, Radial.name))
cbar = plt.colorbar(sc)


plt.show()


#-25.240   -156.750     167.30    -60.320    -35.810     48.530    163.370

# -22.820     70.810     966.30    -50.140     39.040      8.460     90.750