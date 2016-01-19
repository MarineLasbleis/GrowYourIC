#!/usr/local/bin/python
# Time-stamp: <2016-01-19 11:29:07 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Subproject : Getting seismic data (or making them)
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap



def from_seismo_to_cartesian(r, theta, phi):
    """ Calculate the cartesian coordinates from spherical (w/ latitude) coordinates)

    input: 
    r : radius (km)
    theta : latitude (degree)
    phi : longitude (degree)
    output:
    x, y, z
    
    """
    theta = (90-theta) * np.pi/180. # colatitude in rad
    phi = phi * np.pi/180.
    
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z



def from_cartesian_to_seismo(x, y, z):
    """ Calculate the spherical coordinates (w/ latitude) from cartesian coordinates)

    r, theta, phi = from_cartesian_to_seismo(x, y, z)
    input: x, y, z (same length)
    output:
    r : radius (km)
    theta : latitude (degree)
    phi : longitude (degree)
    (same length as the input)

    """

    
    
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arccos(z/r)*180./np.pi # colatitude, in degree

    phi = np.where(y>=0., np.arccos(x/np.sqrt(x**2+y**2)), 2.*np.pi-np.arccos(x/np.sqrt(x**2+y**2)) )
    phi = phi*180./np.pi
    #phi = np.arctan(y/x)*180./np.pi # longitude, in degree
    
    return r, 90.-theta, phi
    


def random_points(type_="turningpoint", seismo="surface"):
    """ Create a random raypath in the IC

    type : type of the output. it can be:
    type="turningpoint" (only the turning point + orientation is given) (default value)
    type="inout" (entrance and exit point of the pathway)

    seismo : distribution in the IC. it can be:
    seismo="surface" : uniform distribution inside an innermost layer of 150km,
    with uniform distribution of orientation (angle zeta)

    """
    depth_min = 0.
    depth_max = 100 #in km
    

    if seismo == "surface":
        depth_turningpoint = np.random.uniform(depth_min, depth_max)
        position = [np.arcsin(np.random.uniform(-1, 1.))*180./np.pi, np.random.uniform(0., 360.)] # theta (latitude) and phi (longitude)
        orientation = np.random.uniform(0., 180.)
    elif seismo == "surface_equatorial":
        depth_turningpoint = np.random.uniform(depth_min, depth_max)
        position = [0., np.random.uniform(0., 360.)] # theta (latitude) and phi (longitude)
        orientation = np.random.uniform(0., 180.)
    else:
        sys.exit("you need to choose a model for the repartition of seismic data")
        
    
    pathway = {'type': type_, 'depth':depth_turningpoint, 'latitude': position[0], 'longitude': position[1], 'angle': orientation}
    return pathway


def raypath_straight(r1, r2, N, coordinate_type="cartesian", Ric=1221):
    """ Calculate the straight line between 2 points
    use cartesian coordinates (change input from spher. to cart. if needed)
    
    r1, r2 : positions of the two points
    N : nombre de points sur la trajectoire
    
    """

    if coordinate_type == "spherical":
        x1, y1, z1 = from_seismo_to_cartesian(Ric-r1['depth'], r1['latitude'], r1['longitude'])
        x2, y2, z2 = from_seismo_to_cartesian(Ric-r2['depth'], r2['latitude'], r2['longitude'])
    elif coordinate_type == "cartesian":
        x1, y1, z1 = r1['x'], r1['y'], r1['z']
        x2, y2, z2 = r2['x'], r2['y'], r2['z']
    else: print "you forget to specify the coordinate type."



    dx = [x2-x1, y2-y1, z2-z1]
    lenght = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2  #between r1 and r2

    trajectory_x = x1 + np.linspace(0, 1, N)*dx[0]
    trajectory_y = y1 + np.linspace(0, 1, N)*dx[1]
    trajectory_z = z1 + np.linspace(0, 1, N)*dx[2]
    dlength = np.sqrt((trajectory_x[1]-trajectory_x[0])**2+ (trajectory_y[1]-trajectory_y[0])**2 +(trajectory_z[1]-trajectory_z[0])**2)


    return trajectory_x, trajectory_y, trajectory_z, dlength, dx
    



def bottom_point_zeta(r, Ric=1221):
    """
    Compute the entry and exit point of a pathway from the bottom point and orientation (angle zeta)
    
    input : position of the bottom point and orientation of the trajectory
    r['type']
    output : entry and exit point
    """
    
    theta = 90 - r['latitude'] # theta is colatitude (degree)
    phi = r['longitude'] # phi is the longitude (degree)
    # perp = [- np.cos(theta), np.sin(theta), 0.]
    radius = Ric-r['depth']
    zeta = r['angle']
    
    if radius < Ric:
        
        dtheta = np.arctan(np.cos(zeta*np.pi/180.)*np.sqrt(Ric**2-radius**2) / radius)*180./np.pi
        theta_1 = theta + dtheta
        theta_2 = theta - dtheta
        dphi = np.arctan(np.sin(zeta*np.pi/180.)*np.sqrt(Ric**2-radius**2) / radius)*180./np.pi
        #np.arctan(np.sin(zeta*np.pi/180.)*np.sqrt(Ric**2-radius**2) / radius)*180./np.pi
        phi_1 = phi + dphi
        phi_2 = phi - dphi

        position_1 = {'type': "spherical", 'depth': 0., 'latitude': 90.-theta_1, 'longitude':phi_1, 'angle': zeta}
        position_2 = {'type': "spherical", 'depth': 0., 'latitude': 90.-theta_2, 'longitude':phi_2, 'angle': zeta}
        
    else:
         print "one of your point is outside the sphere. Ric = ", Ric, " and r = ", radius
         position_1, position_2 = {'type': "spherical", 'depth': 0., 'latitude':theta, 'longitude':phi, 'angle': 0.}
        
    return position_1, position_2



def read_seismic_data(filename, columns, lines=-1):
    """ read seismic data repartition

    input parameters:
    - filename: name of the data file
    - columns: indicates which columns are of interest in the data
    - line: indicates which line will be output (default is -1, all lines are output)
    output:
    - data : columns of the data indicated (number of columns: columns.size, number of lines: line.size)
    """

    data = np.genfromtxt(filename)
    if lines==-1:
        data_ = data[:, columns]
    elif np.size(lines)==1:
        data_ = data[np.ix_([lines], columns)]
    else:
        data_ = data[np.ix_(lines, columns)]
    data = data_
    return data


if __name__ == '__main__':


    N = 50 # number of samples
    RIC = 1221 # inner core radius
    fig, ax = plt.subplots(4)
    fig2, ax2 = plt.subplots(1)
    ax2 = plt.subplot(111, projection='polar')

    for i in range(N):
        pathway = random_points()
        ax[0].scatter(i,pathway['depth'])
        ax[1].scatter(i,pathway['latitude'])
        ax[2].scatter(i,pathway['longitude'])
        ax[3].scatter(i,pathway['angle'])
        ax2.scatter(pathway['longitude']*np.pi/180., RIC-pathway['depth'])


    print from_seismo_to_cartesian(1., 0., 0.), 1, 0, 0
    print from_seismo_to_cartesian(1, 90, 0), 1, 90, 0
    print from_seismo_to_cartesian(1, 0, 90), 1, 0, 90
    print from_seismo_to_cartesian(1, 0, 180), 1, 0, 180
    print from_seismo_to_cartesian(1, -90, 0), 1, -90, 0


    r1 = {'x':0.2, 'y':0.1, 'z':0.8}
    r2 = {'x':1., 'y':1., 'z':1.}
    x, y, z, d, dx = raypath_straight(r1, r2, 10, coordinate_type="cartesian")
    fig3, ax3 = plt.subplots(1)
    ax3.plot(x,y, '-*')
    ax3.plot(y,z, '-+')
    print x, y, z, d, dx

    fig4, ax4 = plt.subplots(1)
    point1 = random_points()
    point2 = random_points()
    point1['latitude'] = 0.
    point2['latitude'] = 0.
    x, y, z, d, dx = raypath_straight(point1, point2, 10, coordinate_type="spherical")
    x1, y1, z1 = from_seismo_to_cartesian(RIC-point1['depth'], point1['latitude'], point1['longitude'])
    x2, y2, z2 = from_seismo_to_cartesian(RIC-point2['depth'], point2['latitude'], point2['longitude'])
    ax4.scatter(x1, y1, c='r')
    ax4.scatter(x2, y2, c='r')
    ax4.plot(x, y)

    # verification raypath with a given point and given angle
    fig5, ax5 = plt.subplots(1,2)
    cm = plt.cm.get_cmap('RdYlBu')
    for i, Angle in enumerate(np.linspace(-90., 90., 10)):
        #(np.linspace(-90., 90., 6)): # Angle is the latitude
        Point =  {'type': "spherical", 'depth': 100., 'latitude':Angle, 'longitude':0., 'angle': 0.}
        pos_1, pos_2 = bottom_point_zeta(Point, Ric=1221)
        ax5[0] = plt.subplot(121, projection='polar')
        ax5[0].scatter(Point['latitude']*np.pi/180., RIC-Point['depth'], marker='x')
        ax5[0].scatter(pos_1['latitude']*np.pi/180., RIC-pos_1['depth'], marker='o')
        ax5[0].scatter(pos_2['latitude']*np.pi/180., RIC-pos_2['depth'], marker='o')
        print Point['latitude'], pos_1['latitude'], pos_2['latitude']
        print Point['longitude'], pos_1['longitude'], pos_2['longitude']
        x, y, z, d, dx = raypath_straight(pos_1, pos_2, 10, coordinate_type="spherical")
        r, theta, phi = from_cartesian_to_seismo(x, y, z)
        ax5[0].scatter(theta*np.pi/180., r, marker='.')
        ax5[1] = plt.subplot(122, projection='polar')
        ax5[1].scatter(Point['longitude']*np.pi/180., RIC-Point['depth'], marker='x')
        ax5[1].scatter(pos_1['longitude']*np.pi/180., RIC-pos_1['depth'], marker='o')
        ax5[1].scatter(pos_2['longitude']*np.pi/180., RIC-pos_2['depth'], marker='o')
        ax5[1].scatter(phi*np.pi/180., r, marker='.')
        print "cartesian"
        print x, y, z
        print "spherical"
        print r, theta, phi


    ## Test Basemap
    
    
    #plt.show()


    fig = plt.subplots()
    m = Basemap(projection='ortho',lat_0=-50,lon_0=30.,resolution='l')
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.25)
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    # draw the edge of the map projection region (the projection limb)
    m.drawmapboundary(fill_color='#99ffff')
    # draw lat/lon grid lines every 30 degrees.
    m.drawmeridians(np.arange(0,360,30))
    m.drawparallels(np.arange(-90,90,30))

    Point =  {'type': "spherical", 'depth': 200., 'latitude': -50., 'longitude':30., 'angle': 0.}
    x, y = m(Point['longitude'], Point['latitude'])
    m.scatter(x, y, 30, marker='o', color='k', zorder=10)
    
    for i, Angle in enumerate(np.linspace(0., 180., 10)):
        Point['angle'] = Angle
        print Point
       
        pos_1, pos_2 = bottom_point_zeta(Point, Ric=1221)
        x1, y1 = m(pos_1['longitude'], pos_1['latitude'])
        x2, y2 = m(pos_2['longitude'], pos_2['latitude'])
        m.scatter(x1, y1, zorder=10)
        m.scatter(x2, y2, zorder=10)
        m.drawgreatcircle(pos_1['longitude'], pos_1['latitude'], pos_2['longitude'], pos_2['latitude'])


    data = read_seismic_data('results.dat', columns=[10, 11, 12, 13], lines=[1,2,3])
    print data
        
    plt.show()
