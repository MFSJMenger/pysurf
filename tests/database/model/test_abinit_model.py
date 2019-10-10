from pytest import fixture
import os
import numpy as np

from pysurf.spp.spp import SurfacePointProvider

@fixture
def input_filename():
    path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(path,'test_abinit_model.inp')

def test_abinit_calc(input_filename):
    if os.path.isfile('db.dat'):
        os.remove('db.dat')
    # set up SPP
    spp = SurfacePointProvider(input_filename)
    coord = spp.refgeo['coord']
    
    # fill DB with three points
    spp.get(spp.refgeo['coord'])
    coord[0, 0] += 1
    spp.get(coord)
    coord[0, 1] += 1
    spp.get(coord)
    coord[0, 2] += 1
    res = spp.get(coord)

    # check interpolation
    coord[0, 2] -= 0.1
    res = spp.get(coord)
    assert(abs(res['energy'][0] - 2.9) < 0.1)
    assert(abs(res['energy'][1] - 0.1) < 0.1)
    refgrad = np.array([[[2, 2, 1.8]], [[ 0. ,  0. , -0.2]]], dtype=float)
    assert(np.amax(np.absolute(res['gradient'] - refgrad)) < 0.1)
