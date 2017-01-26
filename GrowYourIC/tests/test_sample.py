
import py.test 

import GrowYourIC

from GrowYourIC.geodyn_trg import * 
from GrowYourIC.positions import *

#def test_exact_translation_velocity_one():
#    assert exact_translation(SeismoPoint(0,0,0), 1.)==1221.



class TestClass():

    def test_exact_translation_velocity_one(self):
        assert exact_translation(SeismoPoint(0,0,0), 1.)==1221.



class TestAnalytical():
    """ tests the different models in geodynamic.py using analytical solutions"""

    def test_PureGrowth_several_exponents(self):
        for exponent in np.linspace(0.2, 1., 10): 
            Model = PureGrowth()
            Model.set_tauIC(1.)
            Model.set_exponent_growth(exponent)
            Model.set_rICB(1.)
            point = CartesianPoint(0.75, 0., 0.)
            radius = np.linspace(0., 1., 10)
            age = np.zeros_like(radius)

            for ip, r in enumerate(radius):
                r0 = [r, 0., 0.]
                point = positions.CartesianPoint(r, 0., 0.)
                age[ip] = evaluate_singlepoint(point, Model)

            assert np.sqrt(sum((age-(1-(radius)**(1./Model.exponent_growth)))**2.))< 1.e-9
 

