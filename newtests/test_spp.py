import pytest
import numpy as np
from pysurf import SurfacePointProvider


@pytest.fixture
def xtb():
    nstates = 1
    natoms = 3
    return SurfacePointProvider.from_questions(['energy'], nstates, natoms, atomids=['O', 'H', 'H'], config='test.ini')


def test_xtbinterface_energy(xtb):
    res = xtb.request([ [0.00000,       0.00000,       0.11779], 
                  [0.00000,       0.75545,      -0.47116], 
                  [0.00000,      -0.75545,      -0.47116]
                  ], ['energy', 'gradient'])

    assert abs(res['energy'] - -5.07032508030) < 1e-8


def test_xtbinterface_energy_and_gradient(xtb):
    res = xtb.request([ [0.00000,       0.00000,       0.11779], 
                  [0.00000,       0.75545,      -0.47116], 
                  [0.00000,      -0.75545,      -0.47116]
                  ], ['energy', 'gradient'])

    assert abs(res['energy'] - -5.07032508030) < 1e-8
    refgrad = np.array([[2.8722229638859E-17,  -1.0296921891165E-16,   3.6346790530551E-03],
              [-3.5498711199182E-17,  -4.7850233706476E-03,  -1.8173395265276E-03],
              [6.7764815603231E-18,   4.7850233706476E-03,  -1.8173395265275E-03]])
    for e1, r1 in zip(res['gradient'][0].flatten(), refgrad.flatten()):
        assert abs(e1 - r1) < 1e-8
