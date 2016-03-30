#!/usr/local/bin.python
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import pandas as pd

# personal routines
import positions


RICB = 1221.

def read_from_file(filename, names=["station", "PKIKP-PKiKP travel time residual", "zeta", "epicentral distance", "station lat", "stat    ion lon", "event lat", "event lon", "event depth", "in lat", "in lon", "out lat", "out lon", "turn lat", "turn lon", "turn depth", "in    ner core travel time", "PKIKP/PKiKP amplitude ratio"], slices="all"):
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


class SeismicData():
    """ Class for seismic data """
    
    def __init__(self):
        data_points = None
        size = None
        




class FromFile(SeismicData):
    
    def __init__(self, filename="results.dat"):
        
        #seismic data set (from Lauren's file)
        self.filename = filename
        self.slices = ["PKIKP-PKiKP travel time residual", "turn lat", "turn lon", "turn depth", "in lat", "in lon", "out lat", "out lon"]
        self.data = read_from_file(filename, slices=self.slices)
        self.size = self.data.shape

        data_points = []

        for i, row in self.data.iterrows():
            ray = positions.Raypath()
            ray.add_b_t_point(positions.SeismoPoint(RICB-row["turn depth"], row["turn lat"], row["turn lon"]))
            in_Point = positions.SeismoPoint(RICB, row["in lat"], row["in lon"])
            out_Point = positions.SeismoPoint(RICB, row["out lat"], row["out lon"])
            ray.add_in_out(in_Point, out_Point)
            data_points.append(ray)

        
class RandomData(SeismicData):
    pass


if __name__ == '__main__':

    
    f = FromFile()



