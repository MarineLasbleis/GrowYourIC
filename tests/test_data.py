
import py.test 

import GrowYourIC

from GrowYourIC.data import * 

#def test_exact_translation_velocity_one():
#    assert exact_translation(SeismoPoint(0,0,0), 1.)==1221.



class TestClassGeneric():

    def test_SeismicData(self):
        """ assert that empty database is well initialized"""
        data_empty = SeismicData()
        assert data_empty.size==None
        assert data_empty.shortname==None
        assert data_empty.data_points==[]

class TestSeismicDataPoint():

    def test_onepoint(self):
        """ check the construction of a database with 1 point in it """
        data_1pt = RandomData(1)
        assert data_1pt.size == 1
        point = data_1pt.data_points[0].bottom_turning_point

        x, y, z = data_1pt.extract_xyz("bottom_turning_point")
        r, t, p = data_1pt.extract_rtp("bottom_turning_point")
        assert  point.x==x
        assert point.y==y
        assert point.z==z
        assert point.r==r
        assert point.theta==t
        assert point.phi==p


class TestSeismicFromFile():

    def test_load_data(self):
        WD11 = SeismicFromFile('WD11.dat')
        out, err = capfd.readouterr()
        assert out == ""
        

