from pytest import fixture
import os
import numpy as np

from pysurf.spp.spp import SurfacePointProvider

@fixture
def input_filename():
    path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(path,'test.inp')

def test_model(input_filename):
    spp = SurfacePointProvider(input_filename)
    res = spp.get({})
    assert(('energy' in res.keys()) and ('gradient' in res.keys())
           and ('mass' in res.keys()))
    res = spp.get({'coord': np.array([0.0,0.0,0.0]), 
                   'energy': None, 'gradient': None})
    assert(res['energy'].all() == np.array([0.,
                                            0.14479228,
                                            0.17786665]).all())

