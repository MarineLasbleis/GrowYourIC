#!/usr/local/bin/python
# Time-stamp: <2016-03-08 07:40:24 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis


import numpy as np
import pandas as pd

import positions



def read_from_file(filename, names=["station", "PKIKP-PKiKP travel time residual", "zeta", "epicentral distance", "station lat", "station lon", "event lat", "event lon", "event depth", "in lat", "in lon", "out lat", "out lon", "turn lat", "turn lon", "turn depth", "inner core travel time", "PKIKP/PKiKP amplitude ratio"], slices="all"):
    """ read seismic data repartition

    input parameters:
    - filename: name of the data file
    - names: names of the columns for the data set
    - slices: names of columns for the output.
    output:
    - data : pandas DataFrame with all the datas. Columns name are indicated by the variable "names". 
    """
    df = pd.read_table(filename, sep=' ', names=names, skiprows=0)

    if slices != "all":
        df = df[slices]
    return df






if __name__ == '__main__':

    data = read_from_file("results.dat", slices=["turn lat", "turn depth", "turn lon", "in lon"])
    print data.info()

    data_subset = data[["turn lat", "turn depth", "turn lon"]]
    print data_subset.info()
