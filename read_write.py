#!/usr/local/bin/python
# Time-stamp: <2016-02-19 16:16:37 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core

# Author : Marine Lasbleis


import numpy as np

import positions



def read_from_file(filename, columns, lines=-1):
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

    

    print read_from_file("results.dat", columns=[10, 11, 12, 13], lines=[1,2,3])
