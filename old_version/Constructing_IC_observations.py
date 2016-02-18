#!/usr/local/bin/python
# Time-stamp: <2016-01-19 11:28:55 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Subproject : Calculating "synthetics" P-wave velocity from geodynamical models 
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt

import seismo # seismological data repartition
import geodyn # geodynamical models
from mpl_toolkits.basemap import Basemap




def test_translation(N, vt=2., Ric=1221):
    """ Test pure translation

    N : number of points (randomly distributed)
    Ric : radius of IC (km)
    
    """

    fig, ax = plt.subplots(1,3)
    ax[1] = plt.subplot(132, projection='polar')

    for i in range(N):
        pathway = seismo.random_points(seismo="surface")
        ## pathway['depth'] = 0.
        x, y, z = seismo.from_seismo_to_cartesian(Ric-pathway['depth'],
                                                  pathway['longitude'], pathway['longitude'])
        age = geodyn.translation(x/Ric, y/Ric, z/Ric, vt, 1.)
        ax[0].scatter(pathway['longitude'],  age)#, s=age*10)#pathway['depth'])   
        ax[1].scatter(pathway['longitude']*2*np.pi/360, 1221-pathway['depth'], s=age*100)#, s=age*50)
        ax[2].scatter(x, y, s=age*200)
        ## print x/Ric, y/Ric, z/Ric, pathway['depth'], pathway['longitude'], pathway['longitude']
        print age, Ric-pathway['depth'], pathway['longitude'], pathway['longitude']



def test_growth(N, Ric=1221):
    """ Test pure translation

    N : number of points (randomly distributed)
    Ric : radius of IC (km)
    
    """

    fig2, ax = plt.subplots(1)
    fig3, ax2 = plt.subplots(1)
  
    ax2 = plt.subplot(111, projection='polar')
    fig4, ax3 = plt.subplots(1)

    for i in range(N):
        pathway = seismo.random_points(seismo="surface_equatorial")
        ## pathway['depth'] = 0.
        x, y, z = seismo.from_seismo_to_cartesian(Ric-pathway['depth'],
                                                  pathway['longitude'], pathway['longitude'])
        r, age = geodyn.growth(x/Ric, y/Ric, z/Ric, 1.)
        print (1-r)*50
        ax.scatter(pathway['longitude'],  age, s=(1-r)*50)#, s=age*10)#pathway['depth'])   
        ax2.scatter(pathway['longitude']*2*np.pi/360, Ric-pathway['depth'], s=r*10)#, s=age*50)
        ax3.scatter(x, y, s=age*200)
        ## print x/Ric, y/Ric, z/Ric, pathway['depth'], pathway['longitude'], pathway['longitude']

def translation_pathway(N, vt=2., Ric=1221):

    fig, ax = plt.subplots(1)
    fig2, ax2 = plt.subplots(1)
    # use low resolution coastlines.
    map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30))
    map.drawparallels(np.arange(-90,90,30))

    for i in range(N):
        
        Point = seismo.random_points(seismo="surface")
        pos_1, pos_2 = seismo.bottom_point_zeta(Point, Ric=1221)

        x_, y_, z_ = seismo.from_seismo_to_cartesian(Ric-Point['depth'],
                                                  Point['longitude'], Point['latitude'])
        Age = geodyn.translation(x_/Ric, y_/Ric, z_/Ric, vt, 1.)
        age_x, age_y = map(Point['longitude'], Point['latitude'])
        map.scatter(age_x, age_y, s=Age*100, zorder=10)

        #print Age, Ric-Point['depth'], Point['longitude'], Point['longitude']
        Nd = 20
        x, y, z, d, dx = seismo.raypath_straight(pos_1, pos_2, Nd, coordinate_type="spherical")
        r, theta, phi = seismo.from_cartesian_to_seismo(x, y, z)
        
        A = geodyn.translation(x/Ric, y/Ric, z/Ric, vt, 1.)
        Age_average = np.sum(A)/Nd
        
        map.scatter(age_x, age_y, s=Age*100, zorder=20)
        ax.scatter(Point['longitude'], Age_average)
        map.drawgreatcircle(pos_1['longitude'], pos_1['latitude'], pos_2['longitude'], pos_2['latitude'])
        for i, xtheta in enumerate(theta[0:-1]):
            map.drawgreatcircle(phi[i], theta[i], phi[i+1], theta[i+1])
        
        
if __name__ == '__main__':

    N = 10
    
    #test_translation(N)
    # test_growth(N)
    translation_pathway(N)
    plt.show()
    
