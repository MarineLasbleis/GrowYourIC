import py.test 

import GrowYourIC
from GrowYourIC.positions import * 


class TestCoordinateChanges():

    def test_stoc(self):
        assert(from_seismo_to_cartesian(0,0,0)[0]==0.)
        assert(from_seismo_to_cartesian(0,0,0)[1]==0.)
        assert(from_seismo_to_cartesian(0,0,0)[2]==0.)
        assert(from_seismo_to_cartesian(0.,0.,0.)[0]==0.)
        assert(from_seismo_to_cartesian(0.,0.,0.)[1]==0.)
        assert(from_seismo_to_cartesian(0.,0.,0.)[2]==0.)
        np.testing.assert_almost_equal(from_seismo_to_cartesian(1,0,0)[0], 1.)
        np.testing.assert_almost_equal(from_seismo_to_cartesian(1,0,0)[1], 0.)
        np.testing.assert_almost_equal(from_seismo_to_cartesian(1.34,0,0)[0], 1.34)
        np.testing.assert_almost_equal(from_seismo_to_cartesian(1.34,0,0)[2], 0.)
        np.testing.assert_almost_equal(from_seismo_to_cartesian(1,90,0)[2],1.)
        np.testing.assert_almost_equal(from_seismo_to_cartesian(1,0,90)[1],1.)

    def test_CtoS(self):
        np.testing.assert_almost_equal(from_cartesian_to_seismo(np.sqrt(2)/2, np.sqrt(2)/2, 0.)[0],1)
        np.testing.assert_almost_equal(from_cartesian_to_seismo(np.sqrt(2)/2, np.sqrt(2)/2, 0.)[1],0.)
        np.testing.assert_almost_equal(from_cartesian_to_seismo(np.sqrt(2)/2, np.sqrt(2)/2, 0.)[2],45.)
        np.testing.assert_almost_equal(from_cartesian_to_seismo(np.sqrt(2)/2, np.sqrt(2)/2, 1.)[2],45.)

    def test_angular_distance(self):
        np.testing.assert_almost_equal(angular_distance_to_point(0., 0., 0., 0.), 0.)
        np.testing.assert_almost_equal(angular_distance_to_point(12., 13., 12., 13.), 0.)
        np.testing.assert_almost_equal(angular_distance_to_point(0., -100., 0., 0.), 100.)
        np.testing.assert_almost_equal(angular_distance_to_point(-40., 50., 30., 50.), 70.)

    def test_straight_trajectory(self):
        point1 = SeismoPoint(1, 0, 0)
        point2 = SeismoPoint(1, 0, 90)
        trajectory, size = straight_trajectory(point1, point2, 2)
        assert np.array([trajectory]).size==0, "the trajectory should have N-2 points in it (so 0 here)"
        trajectory, size = straight_trajectory(point1, point2, 3)
        assert np.array([trajectory]).size==1, "the trajectory should have N-2 points in it (so 1 here)"
        traj_1 = trajectory[0]
        assert(traj_1.x== 0.5)
        assert(traj_1.y == 0.5)
        np.testing.assert_almost_equal(traj_1.z, 0.)


class TestRaypath(): 

    def test_init(self):
        ray = Raypath()
        assert ray.points == None
        assert ray.direction == None

    def test_add_functions(self):
        ray = Raypath()
        ray.add_b_t_point(SeismoPoint(1.24,0,0))
        ray.add_property({'value':3.})
        assert ray.value== 3.
        ray.add_property({'bottom_turning_point':SeismoPoint(1.24,0,0)})
        assert ray.bottom_turning_point.r == 1.24
        ray.add_property({'bottom_turning_point':SeismoPoint(1,2,3)})
        assert ray.bottom_turning_point.r == 1.24
        ray.add_property({'bottom_turning_point':SeismoPoint(1,2,3)}, brute_force=True)
        assert ray.bottom_turning_point.r == 1.

    def test_straight_in_out1(self):
        ray=Raypath()
        in_point = SeismoPoint(1, 0, 0)
        out_point = SeismoPoint(1, 0, 180)
        ray.add_property({'in_point':in_point, 'out_point':out_point, 'N': 3})
        ray.straight_in_out(ray.N)
        assert ray.length == 2.
        assert np.array([ray.points]).size == 3

    def test_straight_in_out2(self):
        ray=Raypath()
        in_point = SeismoPoint(1, 0, 0)
        out_point = SeismoPoint(0.5, 90, 90)
        ray.add_property({'in_point':in_point, 'out_point':out_point, 'N': 12})
        ray.straight_in_out(ray.N)
        assert ray.length == np.sqrt(1+0.5*0.5)
        assert np.array([ray.points]).size == 12



