from pytest import fixture
import os
import numpy as np

from pysurf.spp.spp import SurfacePointProvider

@fixture
def input_filename():
    path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(path,'test_adc.inp')

def test_abinit_calc(input_filename):
    spp = SurfacePointProvider(input_filename)

    #check properties:
    res = spp.get({})
    assert( ('energy' in res.keys()) and ('gradient' in res.keys()))
    res = spp.get({'coord': spp.refgeo['coord'],
                   'energy': None,
                   'gradient': None})
    #print(res)
    assert np.allclose(res['energy'], np.array([9.99989425e-10, 1.59115397e-01, 1.81781271e-01, 1.91643398e-01]))
