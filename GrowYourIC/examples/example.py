
# # Let's Grow your Own Inner Core!

# ### Choose a model in the list: 
#     - geodyn_trg.TranslationGrowthRotation()
#     - geodyn_static.Hemispheres()
# 
# ### Choose a proxy type:
#     - age
#     - position
#     - phi
#     - theta
#     - growth rate
# 
# ### set the parameters for the model : geodynModel.set_parameters(parameters)
# ### set the units : geodynModel.define_units()
# 
# ### Choose a data set:
#     - data.SeismicFromFile(filename) # Lauren's data set
#     - data.RandomData(numbers_of_points)
#     - data.PerfectSamplingEquator(numbers_of_points)
#         organized on a cartesian grid. numbers_of_points is the number of points along the x or y axis. The total number of points is numbers_of_points**2*pi/4
#         - as a special plot function to show streamlines: plot_c_vec(self,modelgeodyn)
#     - data.PerfectSamplingEquatorRadial(Nr, Ntheta)
#         same than below, but organized on a polar grid, not a cartesian grid.
# 
# 
# ### Extract the info:
#     - calculate the proxy value for all points of the data set: geodyn.evaluate_proxy(data_set, geodynModel)
#     - extract the positions as numpy arrays: extract_rtp or extract_xyz
#     - calculate other variables: positions.angular_distance_to_point(t,p, t_point, p_point)


# import statements
import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math

from GrowYourIC import positions, geodyn, geodyn_trg, geodyn_static, plot_data, data

plt.rcParams['figure.figsize'] = (8.0, 3.0) #size of figures
cm = plt.cm.get_cmap('viridis')
cm2 = plt.cm.get_cmap('winter')


# ## Define the geodynamical model

geodynModel = geodyn_trg.TranslationGrowthRotation() #can do all the models presented in the paper
# geodynModel = geodyn_static.Hemispheres() #this is a static model, only hemispheres. 


# Change the values of the parameters to get the model you want (here, parameters for .TranslationGrowthRotation())


age_ic_dim = 1e9 #in years
rICB_dim = 1221. #in km
translation_velocity_dim = 4e-10 #m.s, value for today's Earth with Q_cmb = 10TW (see Alboussiere et al. 2010)
time_translation = rICB_dim*1e3/translation_velocity_dim /(np.pi*1e7)
maxAge = 2.*time_translation/1e6
print("The translation recycles the inner core material in {0:.2f} million years".format(maxAge))
print("Translation velocity is {0:.2e} km/years".format(translation_velocity_dim*np.pi*1e7/1e3))

units = None #we give them already dimensionless parameters. 
rICB = 1.
age_ic = 1.
omega = 0. #-0.5*np.pi # Rotation rates has to be in ]-np.pi, np.pi[
velocity_amplitude = translation_velocity_dim*age_ic_dim*np.pi*1e7/rICB_dim/1e3
velocity_center = [0., 100.]#center of the eastern hemisphere
velocity = geodyn_trg.translation_velocity(velocity_center, velocity_amplitude)
exponent_growth = 0.3 


# Define a proxy type, and a proxy name (to be used in the figures to annotate the axes)
# 
# You can re-define it later if you want (or define another proxy_type2 if needed)


proxy_type = "age"
proxy_name = "age (Myears)"
proxy_lim = [0, maxAge] #or None

fig_name = "figures/test_" #to name the figures


# ### Parameters for the geodynamical model
# 
# This will input the different parameters in the model.

parameters = {'units': units,
              'rICB': rICB, 
              'tau_ic':age_ic,
              'vt': velocity,
              'exponent_growth': 0.5,
              'omega': omega,
              'proxy_type': proxy_type}
geodynModel.set_parameters(parameters)
geodynModel.define_units()


# ## Different data set and visualisations

# ### Perfect sampling at the equator (to visualise the flow lines)
# 
# You can add more points to get a better precision.

npoints = 50 #number of points in the x direction for the data set. 
data_set = data.PerfectSamplingEquator(npoints, rICB = 1.)
data_set.method = "bt_point"
proxy = geodyn.evaluate_proxy(data_set, geodynModel, proxy_type="age", verbose = False)
data_set.plot_c_vec(geodynModel, proxy=proxy, cm=cm, nameproxy="age (Myears)")
plt.savefig(fig_name+"equatorial_plot.pdf", bbox_inches='tight')


# ### Perfect sampling in the first 100km (to visualise the depth evolution)

data_meshgrid = data.Equator_upperpart(150,150)
data_meshgrid.method = "bt_point"
proxy_meshgrid = geodyn.evaluate_proxy(data_meshgrid, geodynModel, proxy_type=proxy_type, verbose = False)
r, t, p = data_meshgrid.extract_rtp("bottom_turning_point")

fig3, ax3 = plt.subplots(figsize=(8, 2))
X, Y, Z  = data_meshgrid.mesh_RPProxy(proxy_meshgrid)
sc = ax3.contourf(Y, rICB_dim*(1.-X), Z, 100, cmap=cm)
sc2 = ax3.contour(sc, levels=sc.levels[::15], colors = "k")
ax3.set_ylim(-0, 120)
fig3.gca().invert_yaxis()
ax3.set_xlim(-180,180)
cbar = fig3.colorbar(sc)
#cbar.set_clim(0, maxAge)
cbar.set_label(proxy_name)
ax3.set_xlabel("longitude")
ax3.set_ylabel("depth below ICB (km)")

plt.savefig(fig_name+"meshgrid.pdf", bbox_inches='tight')


# ### Random data set, in the first 100km  - bottom turning point only
# 
# #### Calculate the data

# random data set
data_set_random = data.RandomData(3000)
data_set_random.method = "bt_point"

proxy_random = geodyn.evaluate_proxy(data_set_random, geodynModel, proxy_type=proxy_type, verbose=False)
#if proxy_type == "age":
    ## domain size and Vp
#    proxy_random_size = geodyn.evaluate_proxy(data_set_random, geodynModel, proxy_type="domain_size", verbose=False)
#    proxy_random_dV = geodyn.evaluate_proxy(data_set_random, geodynModel, proxy_type="dV_V", verbose=False)


r, t, p = data_set_random.extract_rtp("bottom_turning_point")
dist = positions.angular_distance_to_point(t, p, *velocity_center)

## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy_random,s=8, zorder=10, cmap=cm, edgecolors='none')
plt.title("Dataset: {},\n geodynamic model: {}".format(data_set_random.name, geodynModel.name))
cbar = plt.colorbar(sc)
cbar.set_label(proxy_name)
#fig.savefig(fig_name+data_set_random.name+"map.pdf", bbox_inches='tight')

## phi and distance plots
fig, ax = plt.subplots(1,2, figsize=(8.0, 3.0))
sc1 = ax[0].scatter(p, proxy_random, c=abs(t),s=3, cmap=cm2, vmin =-0, vmax =90, linewidth=0)
phi = np.linspace(-180,180, 50)
analytic_equator = np.maximum(2*np.sin((phi-10)*np.pi/180.)*rICB_dim*1e3/translation_velocity_dim /(np.pi*1e7)/1e6,0.)
ax[0].plot(phi,analytic_equator, 'r', linewidth=2)
ax[0].set_xlabel("longitude")
ax[0].set_ylabel(proxy_name)
if proxy_lim is not None:
    ax[0].set_ylim(proxy_lim)
sc2 = ax[1].scatter(dist, proxy_random, c=abs(t), cmap=cm2, vmin=-0, vmax =90, s=3, linewidth=0)
ax[1].set_xlabel("angular distance to ({}, {})".format(*velocity_center))
phi = np.linspace(-90,90, 100)
if proxy_type == "age":
    analytic_equator = np.maximum(2*np.sin((-phi)*np.pi/180.)*rICB_dim*1e3/translation_velocity_dim /(np.pi*1e7)/1e6,0.)
    ax[1].plot(phi+90,analytic_equator, 'r', linewidth=2)
ax[1].set_xlim([0,180])
ax[0].set_xlim([-180,180])
cbar = fig.colorbar(sc1)
cbar.set_label("longitude: abs(theta)")
if proxy_lim is not None:
    ax[1].set_ylim(proxy_lim)
## figure with domain size and Vp

#fig.savefig(fig_name +data_set_random.name+ 'long_dist.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(8, 2))
sc=ax.scatter(p,rICB_dim*(1.-r), c=proxy_random, s=10,cmap=cm, linewidth=0)
ax.set_ylim(-0,120)
fig.gca().invert_yaxis()
ax.set_xlim(-180,180)
cbar = fig.colorbar(sc)
if proxy_lim is not None:
    cbar.set_clim(0, maxAge)
ax.set_xlabel("longitude")
ax.set_ylabel("depth below ICB (km)")
cbar.set_label(proxy_name)

#fig.savefig(fig_name+data_set_random.name+"depth.pdf", bbox_inches='tight')


# ### Real Data set from Waszek paper

## real data set
data_set = data.SeismicFromFile("../GrowYourIC/data/WD11.dat")
data_set.method = "bt_point"
geodynModel.proxy_type = "age"
proxy2 = geodyn.evaluate_proxy(data_set, geodynModel, proxy_type="age", verbose=False)


r, t, p = data_set.extract_rtp("bottom_turning_point")
dist = positions.angular_distance_to_point(t, p, *velocity_center)

## map
m, fig = plot_data.setting_map() 
x, y = m(p, t)
sc = m.scatter(x, y, c=proxy2,s=8, zorder=10, cmap=cm, edgecolors='none')
plt.title("Dataset: {},\n geodynamic model: {}".format(data_set_random.name, geodynModel.name))
cbar = plt.colorbar(sc)
cbar.set_label(proxy_name)
#fig.savefig(fig_name+data_set_random.name+"map.pdf", bbox_inches='tight')

## phi and distance plots
fig, ax = plt.subplots(1,2, figsize=(8.0, 3.0))
sc1 = ax[0].scatter(p, proxy2, c=abs(t),s=3, cmap=cm2, vmin =-0, vmax =90, linewidth=0)
phi = np.linspace(-180,180, 50)
analytic_equator = np.maximum(2*np.sin((phi-10)*np.pi/180.)*rICB_dim*1e3/translation_velocity_dim /(np.pi*1e7)/1e6,0.)
ax[0].plot(phi,analytic_equator, 'r', linewidth=2)
ax[0].set_xlabel("longitude")
ax[0].set_ylabel(proxy_name)
if proxy_lim is not None:
    ax[0].set_ylim(proxy_lim)
sc2 = ax[1].scatter(dist, proxy2, c=abs(t), cmap=cm2, vmin=-0, vmax =90, s=3, linewidth=0)
ax[1].set_xlabel("angular distance to ({}, {})".format(*velocity_center))
phi = np.linspace(-90,90, 100)
if proxy_type == "age":
    analytic_equator = np.maximum(2*np.sin((-phi)*np.pi/180.)*rICB_dim*1e3/translation_velocity_dim /(np.pi*1e7)/1e6,0.)
    ax[1].plot(phi+90,analytic_equator, 'r', linewidth=2)
ax[1].set_xlim([0,180])
ax[0].set_xlim([-180,180])
cbar = fig.colorbar(sc1)
cbar.set_label("longitude: abs(theta)")
if proxy_lim is not None:
    ax[1].set_ylim(proxy_lim)
## figure with domain size and Vp

#fig.savefig(fig_name +data_set_random.name+ 'long_dist.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(8, 2))
sc=ax.scatter(p,rICB_dim*(1.-r), c=proxy2, s=10,cmap=cm, linewidth=0)
ax.set_ylim(-0,120)
fig.gca().invert_yaxis()
ax.set_xlim(-180,180)
cbar = fig.colorbar(sc)
if proxy_lim is not None:
    cbar.set_clim(0, maxAge)
ax.set_xlabel("longitude")
ax.set_ylabel("depth below ICB (km)")
cbar.set_label(proxy_name)


plt.show()
