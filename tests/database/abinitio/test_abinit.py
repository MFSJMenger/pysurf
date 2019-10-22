from pytest import fixture
import os
import numpy as np

from pysurf.spp.spp import SurfacePointProvider

@fixture
def input_filename():
    path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(path,'test_abinit.inp')

def test_abinit_calc(input_filename):
    spp = SurfacePointProvider(input_filename)
    res = spp.get({})
    assert(('energy' in res.keys()) and ('gradient' in res.keys()))
    
    res = spp.get({'coord': spp.refgeo['coord'],
                   'energy': None, 'gradient': None})

    assert(res['energy'].all() == np.array([-8.00014277e-10,
                                            -1.52319380e-02,
                                            -1.40437478e-02, 
                                            1.29286264e-02]).all())
