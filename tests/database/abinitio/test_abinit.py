from pytest import fixture
import os

from pysurf.spp.spp import SurfacePointProvider

@fixture
def input_filename():
    path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(path,'test_abinit.inp')

def test_abinit_calc(input_filename):
    spp = SurfacePointProvider(input_filename)
    print(spp.db.get('coord',1))
