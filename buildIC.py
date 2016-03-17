#!/usr/local/bin/python
# Time-stamp: <2016-03-08 07:42:35 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps

#personal routines
import positions
import plot_data
import read_write
import geodynamic



def Translation_analytic_solution():

    age_IC = 1.e6 #in years
    velocity_translation = 2.01*1221./1.e6 #in km/years 

    #Initialize map
    m, fig = plot_data.setting_map()
    cm = plt.cm.get_cmap('RdYlBu')
    
    #seismic data set (from Lauren's file)
    data_points = read_write.read_from_file("results.dat", columns=[1, 13, 14, 15], lines=-1)
    nlines, ncolumns = data_points.shape

    #translate it in the correct format
    dataset = []
    translation_dataset = []
    phi, theta, age = [], [], []
    print nlines, "points to write."
    for i in range(nlines):
        if i%100==0: #print every 100 values
            print "Writing point ",  i, ". Coordinates: ", data_points[i, :]
        # translate data in the correct format
        dataset.append(positions.Raypath("BT-point",
                                         1221.-data_points[i,3],
                                         data_points[i, 1],
                                         data_points[i, 2]))
        # set up the geodynamical model
        # because it's analytical solution, we don't need to define a grid
        translation_dataset.append(geodynamic.Translation(positions.Point(1221.-data_points[i,3],
                                                          data_points[i, 1],
                                                          data_points[i, 2],
                                                          "seismo")))
        translation_dataset[i].analytical(velocity_translation, age_IC)


        phi.append(translation_dataset[i].initial_position.phi)
        theta.append(translation_dataset[i].initial_position.theta)
        age.append(translation_dataset[i].exact_solution)

    
    x, y = m(phi, theta)
    m.scatter(x, y, c=age, zorder=10, cmap=plt.cm.RdYlGn)
    m.colorbar()
    
    plt.show()
    
    ## x, y = m(phi, theta)
    ## m.scatter(x, y, c=(data_points[:,0]-min(data_points[:,0]))*1, zorder=10, cmap=plt.cm.RdYlGn)
    ## m.colorbar()

    ## plt.show()



if __name__ == '__main__':


    Translation_analytic_solution()


    # build IC model
    RIC = 1221.
    age_IC = 1.e6 #in years
    velocity_translation = 2.01*1221./1.e6 #in km/years 

    #Initialize map
    m, fig = plot_data.setting_map()
    cm = plt.cm.get_cmap('RdYlBu')


    # points coverage: only surface now
    Ntheta, Nphi = 100, 100
    Ntotal = Ntheta*Nphi
    theta = np.linspace(-90, 90, Ntheta)
    phi = np.linspace(0, 360, Nphi)

    # Geodynamical model definition
    model = []
    for it in range(Ntheta):
        for ip in range(Nphi):
            model.append(geodynamic.Translation(positions.Point(RIC,
                                                                theta[it],
                                                                phi[ip],
                                                                "seismo")))
    
