#!/usr/local/bin/python
# Time-stamp: <2016-01-19 11:29:03 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Subproject : Geodynamic models
# Author : Marine Lasbleis




import numpy as np
import matplotlib.pyplot as plt


def translation(x, y, z, vt, t0=1.):
    """ Pure translation, constant (positive) velocity along x-axis

    the analytical solution for age inside a sphere is
    t = (-(1-y^2-z^2)^0.5-x)/vt +t0
    t0 is the (global) age of the inner core,
    x,y,z the coordinates of the interesting point. 
    """

    return (-(1-y**2-z**2)**0.5-x)/vt +t0




def growth(x, y, z, t0, alpha = 0.5):
    """ Pure growth for R~t^alpha

    the analytical solution for this is
    t = (t0**alpha + ln(r))**(1/alpha)
    t0 is the (global) age of the inner core,
    r is the coordinates of the interesting point
    r = sqrt(x**2+y**2+z**2)
    """
    r = np.sqrt(x**2+y**2+z**2)
    
    return r, (t0**alpha + np.log(r))**(1/alpha)
    

if __name__ == '__main__':


    fig, ax = plt.subplots(1,2)
    
    x = np.linspace(-1, 1, 100)
    z = np.linspace(-1, 1, 100)

    X, Z = np.meshgrid(x, z)
    R = X**2. + Z**2
    AgeTranslation = np.zeros(R.shape)
    AgeGrowth = np.zeros(R.shape)

    
    for i, radius in enumerate(R.flat):
        if radius < 1 :
            AgeTranslation.flat[i] = translation(X.flat[i], 0.,  Z.flat[i], 0.1, 1)
            AgeGrowth.flat[i] = growth(X.flat[i], 0.,  Z.flat[i], 1)

    ax[0].contourf(X, Z, AgeTranslation, 10)
    ax[0].axis('equal')
    ax[1].contourf(X, Z, AgeGrowth, 10)
    ax[1].axis('equal')

    x = np.linspace(0.01, 1, 20)
    fig2, ax2 = plt.subplots(1)
    ax2.plot(x, growth(x, 0., 0., 1))
    
    plt.show()
