#!/usr/local/bin/python
# Time-stamp: <2016-01-19 11:45:37 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Subproject : Constructing IC observations
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import yaml

import seismo # seismological data repartition
import geodyn # geodynamical models

def seismic_data(Method, N=20):
    """ Give the position of the point and raypath of interest,
    given a method (random sampling, existing dataset, etc.)

    """

    if Method['s_type']=='random':
        Point = seismo.random_points(seismo=Method['random_type'])
        In, Out = seismo.bottom_point_zeta(Point, Ric=1221)

    elif Method['s_type']=='real data':
        data = seismo.read_seismic_data(Method['fname'], Method['columns'], lines=Method['lines'])
        Point = {'type': "spherical", 'depth':data[:,6], 'latitude':data[:,4], 'longitude':data[:,5], 'angle': 0.}
        In = {'type': "spherical", 'depth':0., 'latitude':data[:,0], 'longitude':data[:,1], 'angle': 0.}
        Out = {'type': "spherical", 'depth':0., 'latitude':data[:,2], 'longitude':data[:,3], 'angle': 0.}
    
    Points = {'turning point':Point, 'in':In, 'out':Out}
    x, y, z, d, dx = seismo.raypath_straight(In, Out, N, coordinate_type="spherical")
    r, theta, phi = seismo.from_cartesian_to_seismo(x, y, z)
    Trajectory = {'type': "spherical", 'depth':Ric-r, 'latitude':theta, 'longitude':phi, 'angle': Point['angle']}

    return Points, Trajectory



def geodynamics(Points, Trajectory, Method, Ric=1221):
    """ Compute the value of the proxy over the trajectory
    
    """

    Nd = Trajectory['latitude'].size 
    
    if Method['g_type']=='translation':

        velocity = Method['v']
        Point = Points['turning point']
        x_, y_, z_ = seismo.from_seismo_to_cartesian(Ric-Point['depth'],
                                                Point['latitude'], Point['longitude'])
        x, y, z = seismo.from_seismo_to_cartesian(Ric-Trajectory['depth'],
                                                Trajectory['latitude'], Trajectory['longitude'])
        Age = geodyn.translation(x_/Ric, y_/Ric, z_/Ric, velocity, 1.)
        A = geodyn.translation(x/Ric, y/Ric, z/Ric, velocity, 1.)
        Age_average = np.sum(A)/Nd

    return Age_average, Age



if __name__ == '__main__':


    with open("param.dat", 'r') as stream:
        PARAM = yaml.load(stream)
    g_method = PARAM['geodynamics']
    s_method = PARAM['seismic_method']
    Ric = PARAM['Ric']
    Nd = PARAM['Nd']
    
    ## vt=3. # translation velocity
    ## Ric=1221. # radius IC
    ## Nd = 20 #number of points in the trajectory

    ## N = 3184 #number of points 3184
    
    ## g_method = {'g_type':'translation', 'v':vt}
    ## #s_method = {'s_type':'random', 's_N':10, 'random_type':'surface'}
    ## s_method = {'s_type':'real data', 'fname':'results.dat', 'columns':[9, 10, 11, 12, 13, 14, 15], 'lines':1}
    g_method.update(s_method)
    method = g_method
    
    # Initialisation of figures
    fig, ax = plt.subplots(1)
    fig2, ax2 = plt.subplots(1)
    # use low resolution coastlines.
    map_data = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    # draw coastlines, country boundaries, fill continents.
    map_data.drawcoastlines(linewidth=0.25)
    map_data.drawcountries(linewidth=0.25)
    map_data.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map_data.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.
    map_data.drawmeridians(np.arange(0,360,30))
    map_data.drawparallels(np.arange(-90,90,30))

    
    
    # Choose a point + compute raypath (or obtain raypath from data set)
   
    
    for line in range(s_method['Nlines']):

        method['lines'] = line

        Points, Trajectory = seismic_data(method, N=Nd)

        # Choose a model and calculate the proxy (average over the raypath)
        Age_average, Age = geodynamics(Points, Trajectory, method)

        # Plot the calculated points
        age_x, age_y = map_data(Points['turning point']['longitude'], Points['turning point']['latitude'])
        map_data.scatter(age_x, age_y, s=Age*100, zorder=10)
        map_data.scatter(age_x, age_y, s=Age*100, zorder=20)
        map_data.drawgreatcircle(Points['in']['longitude'], Points['in']['latitude'], Points['out']['longitude'], Points['out']['latitude'])
        ax.scatter(Points['turning point']['longitude'], Age_average)

    
    # Save the data


    # Print the plot
    plt.show()
