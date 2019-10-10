from pytest import fixture
import os

from pysurf.spp.spp import SurfacePointProvider

@fixture
def input_filename():
    path = os.path.dirname(os.path.realpath(__file__))
    return os.path.join(path,'test_abinit_interpolation.inp')

def test_abinit_calc(input_filename):
    spp = SurfacePointProvider(input_filename)
#    print(spp.get(spp.refgeo['coord']))
#    print(spp.db.get('coord',0))
#    print(spp.db.get('gradient',0))
    spp.get(spp.refgeo['coord'])
    coord = spp.refgeo['coord']
    coord[0] += 0.2
    print(coord)
    spp.get(coord)
    coord[1] += 0.3
    print(coord)
    spp.get(coord)
    spp.get(coord)
    print(spp.dbinter.db.get_dimension_size('frame'))
#    assert(spp.db.get_dimension_size('frames') == 1)
