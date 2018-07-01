from __future__ import absolute_import

# import statements
import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math

import GrowYourIC

from GrowYourIC import positions
from GrowYourIC import geodyn, geodyn_trg, geodyn_static
from GrowYourIC import plot_data
from GrowYourIC import data

plt.rcParams['figure.figsize'] = (8.0, 3.0) #size of figures
cm = plt.cm.get_cmap('summer')
cm2 = plt.cm.get_cmap('winter')


geodynModel = geodyn_trg.TranslationGrowthRotation()

age_ic_dim = 1e9 #in years
rICB_dim = 1221. #in km
velocity_center = [0., 100.]#center of the eastern hemisphere
center = [0,-80] #center of the western hemisphere
units = None #we give them already dimensionless parameters. 
rICB = 1.
age_ic = 1.

v_fast = 10.3
v_dim = 12.6
omega_fast = 7.85
time_translation = rICB_dim*1e3/4e-10/(np.pi*1e7)
maxAge = 2.*time_translation/1e6
velocity_fast = geodyn_trg.translation_velocity(velocity_center, v_fast)
exponent_fast = 0.1

proxy_type = "age"
proxy_name = "age (Myears)" #growth rate (km/Myears)"
proxy_lim = [0, maxAge]

print("The translation recycles the inner core material in {0:.2f} million years".format(maxAge))

parameters = dict({'units': units,
              'rICB': rICB, 
              'tau_ic':age_ic,
              'vt': velocity_fast,
              'exponent_growth': exponent_fast,
              'omega': 0.,
              'proxy_type': proxy_type, 
              'proxy_name': proxy_name,
              'proxy_lim': proxy_lim})
geodynModel.set_parameters(parameters)
geodynModel.define_units()



## Visualize the flow (equatorial cross section)

npoints = 30 #number of points in the x direction for the data set. 
data_set = data.PerfectSamplingEquator(npoints, rICB = 1.)
data_set.method = "bt_point"
proxy_ = geodyn.evaluate_proxy(data_set, geodynModel, proxy_type="age", verbose = False)
data_set.plot_c_vec(geodynModel, proxy=proxy_, cm=cm, nameproxy="age (Myears)")


npoints = 30 #number of points in the x direction for the data set. 
data_set = data.PerfectSamplingSurface(npoints, rICB = 1., depth=0.01)
data_set.method = "bt_point"
surface1 = geodyn.evaluate_proxy(data_set, geodynModel, verbose = False)
X, Y, Z  = data_set.mesh_TPProxy(surface1)
m, fig = plot_data.setting_map()
y, x = m(Y, X)
sc = m.contourf(y, x, Z, 30, cmap=cm, zorder=2, edgecolors='none')
cbar = plt.colorbar(sc)
cbar.set_label(geodynModel.proxy_name)


# perfect repartition in depth (for meshgrid plots)
data_meshgrid = data.Equator_upperpart(30,30)
data_meshgrid.method = "bt_point"
meshgrid1 = geodyn.evaluate_proxy(data_meshgrid, geodynModel, verbose = False)
fig3, ax3 = plt.subplots(figsize=(8, 2))
X, Y, Z  = data_meshgrid.mesh_RPProxy(meshgrid1)
sc = ax3.contourf(Y, rICB_dim*(1.-X), Z, 100, cmap=cm)
sc2 = ax3.contour(sc, levels=sc.levels[::15], colors = "k")
ax3.set_ylim(-0, 120)
fig3.gca().invert_yaxis()
ax3.set_xlim(-180,180)
cbar = fig3.colorbar(sc)
#cbar.set_clim(0, maxAge)
cbar.set_label(geodynModel.proxy_name)
ax3.set_xlabel("longitude")
ax3.set_ylabel("depth below ICB (km)")



## real data set - WD13
data_set = data.SeismicFromFile("../GrowYourIC/data/WD11.dat")
data_set.method = "bt_point"
proxy1 = geodyn.evaluate_proxy(data_set, geodynModel, verbose=False)

r, t, p = data_set.extract_rtp("bottom_turning_point")
dist = positions.angular_distance_to_point(t, p, *center)

## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy1,s=8, zorder=10, cmap=cm, edgecolors='none')
cbar = plt.colorbar(sc)
cbar.set_label(geodynModel.proxy_name)

fig, ax = plt.subplots(figsize=(8, 2))
sc=ax.scatter(p,rICB_dim*(1.-r), c=proxy1, s=10,cmap=cm, linewidth=0)
ax.set_ylim(-0,120)
fig.gca().invert_yaxis()
ax.set_xlim(-180,180)
cbar = fig.colorbar(sc)
if proxy_lim is not None:
    cbar.set_clim(0, maxAge)
ax.set_xlabel("longitude")
ax.set_ylabel("depth below ICB (km)")
cbar.set_label(geodynModel.proxy_name)

## phi and distance plots
fig, ax = plt.subplots(1,1, figsize=(4.0, 2.5))
sc1 = ax.scatter(p, proxy1, c=abs(t),s=3, cmap=cm2, vmin =-0, vmax =90, linewidth=0)
phi = np.linspace(-180,180, 50)
ax.set_xlabel("longitude")
ax.set_ylabel(proxy_name)
if proxy_lim is not None:
    ax.set_ylim(proxy_lim)
phi = np.linspace(-90,90, 100)
if proxy_type == "age":
    analytic_equator = np.maximum(2*np.sin((phi-10)*np.pi/180.)*rICB_dim*1e3/v_dim*1e-3,0.)
    ax.plot(phi,analytic_equator, 'r', linewidth=2)
ax.set_xlim([-180,180])
cbar = fig.colorbar(sc1)
cbar.set_label("longitude: abs(theta)")



# random data set -
data_set_random = data.RandomData(100)
data_set_random.method = "bt_point"
proxy_random1 = geodyn.evaluate_proxy(data_set_random, geodynModel, verbose=False)

r, t, p = data_set_random.extract_rtp("bottom_turning_point")
dist = positions.angular_distance_to_point(t, p, *center)

## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy_random1,s=8, zorder=10, cmap=cm, edgecolors='none')
cbar = plt.colorbar(sc)
cbar.set_label(geodynModel.proxy_name)

fig, ax = plt.subplots(figsize=(8, 2))
sc=ax.scatter(p,rICB_dim*(1.-r), c=proxy_random1, s=10,cmap=cm, linewidth=0)
ax.set_ylim(-0,120)
fig.gca().invert_yaxis()
ax.set_xlim(-180,180)
cbar = fig.colorbar(sc)
if proxy_lim is not None:
    cbar.set_clim(0, maxAge)
ax.set_xlabel("longitude")
ax.set_ylabel("depth below ICB (km)")
cbar.set_label(geodynModel.proxy_name)

## phi and distance plots
fig, ax = plt.subplots(1,1, figsize=(4.0, 2.5))
sc1 = ax.scatter(p, proxy_random1, c=abs(t),s=3, cmap=cm2, vmin =-0, vmax =90, linewidth=0)
phi = np.linspace(-180,180, 50)
ax.set_xlabel("longitude")
ax.set_ylabel(proxy_name)
if proxy_lim is not None:
    ax.set_ylim(proxy_lim)
phi = np.linspace(-90,90, 100)
if proxy_type == "age":
    analytic_equator = np.maximum(2*np.sin((phi-10)*np.pi/180.)*rICB_dim*1e3/v_dim*1e-3,0.)
    ax.plot(phi,analytic_equator, 'r', linewidth=2)
ax.set_xlim([-180,180])
cbar = fig.colorbar(sc1)
cbar.set_label("longitude: abs(theta)")



plt.show()
