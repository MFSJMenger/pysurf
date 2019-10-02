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
    print(spp.get(np.array([0.0,0.0,0.0])))
