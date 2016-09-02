#!/usr/local/bin/python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis

# Seismic properties of material as function of ka (adimensional frequency): P and S wave velocity and attenuation.
# Please refer to Calvet and Margerin 2008 (figures 3 and 4)

from __future__ import division
from __future__ import absolute_import


import numpy as np
import matplotlib.pyplot as plt  # for figures
import warnings
import scipy.io as sio

# IMPORTANT : all the functions here use dimensional quantities in unit SI !


def export_matlab_data(name_data, file_name="./CM2008/data.mat"):
    """ Values for polynomial fit of the velocity and attenuation are included in data.mat. Scipy export them as a dictionary."""
    return sio.loadmat(file_name)[name_data]


def domain_size(age):
    """ domain size of a crystal in the inner core.

    age is dimensionnal (in seconds)
    Value of M \cdot \gamma = 8.10^{10} m2/s from Bergman 2010 (see Geballe 2013, paragraph 12)
    output is the domain size, given in meters.
    """
    Mgamma = 8.e-10  # m^2/s grain boundary mobility * surface energy
    # exponent: 0.5 for pure materials, 1/3 to 1/4 for alloys
    return np.sqrt(Mgamma * age)


def adimensional_frequency(size, v=11030., freq=1.):
    """ k_0 a = domain_size * 2 pi freq / v """
    return size * 2. * np.pi * freq / v


def convert_CM2008_velocity(kR, poly):
    """ Data from Fig 3 and 4 of Calvet and Margerin 2008 have been stored as polynomial values. This function transform them to a function:
        for example, with poly = Vp_poly, 
        give Vp as function of kR (adimensional frequency)"""
#    if any(min(abs(np.log10(kR)+2),abs(np.log10(kR-1))) < 1e-4):
#        warnings.warn('heaviside fn may be computing 0.5')
    velocity = heaviside(kR - 10.**(-2.)) * heaviside(10. - kR) * np.polyval(poly, np.log10(kR)) +\
        heaviside(10.**(-2) - kR) * np.polyval(poly, -2.) +\
        heaviside(kR - 10.) * np.polyval(poly, 1.)
    return velocity


def convert_CM2008_attenuation(kR, poly):
    """ Data from Fig 3 and 4 of Calvet and Margerin 2008 have been stored as polynomial values. This function transform them to a function:
        for example, with poly = Qp_poly, 
        give attenuation as function of kR (adimensional frequency)"""
#    if min(abs(np.log10(kR)+2),abs(np.log10(kR-1))) < 1e-4:
#        warnings.warn('heaviside fn may be computing 0.5')
    velocity = heaviside(kR - 10**-1) * heaviside(10**1 - kR) * np.polyval(poly, np.log10(kR)) +\
        heaviside(10**-1 - kR) * np.polyval(poly, -1) +\
        heaviside(kR - 10**1) * np.polyval(poly, 1)
    velocity = 10**velocity
    return velocity


def heaviside(x):
    """ numpy does not define an heaviside function.Let's define it here."""
    return 0.5 * (np.sign(x) + 1)


if __name__ == "__main__":

    x = 10**(np.linspace(-4, 4, 200))
    plt.plot(np.log10(x), convert_CM2008_velocity(
        x, export_matlab_data("Belanoshko_Vp_poly")), 'o-')

    plt.show()
