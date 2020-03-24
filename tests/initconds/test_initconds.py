from pytest import fixture
import os
import numpy as np



from pysurf.sampling.initialconditions import InitialConditionsFactory
from pysurf.utils.osutils import exists_and_isfile


pyrazine_crd = np.array([[ 2.13131002e+00,  1.31728882e+00,  9.30480000e-04],
                         [6.90960000e-04,  2.66329754e+00,  2.31798000e-03],
                         [-2.13078722e+00,  1.31820446e+00,  1.22076000e-03],
                         [-2.13142558e+00, -1.31717271e+00, -8.40390000e-04],
                         [-3.90820150e+00,  2.37650649e+00,  2.09824000e-03],
                         [-3.90935032e+00, -2.37461841e+00, -1.87133000e-03],
                         [-6.04090000e-04, -2.66329813e+00, -2.43821000e-03],
                         [  2.13066752e+00, -1.31832153e+00, -1.20462000e-03],
                         [ 3.90824960e+00, -2.37637730e+00, -2.09090000e-03],
                         [ 3.90941282e+00,  2.37447189e+00,  1.93468000e-03]])

@fixture
def filepath_molden():
    return 'test_molden.inp'

@fixture
def filepath_freq():
    return 'test_freq.inp'

def delete_database(filename):
    if exists_and_isfile(filename):
        os.remove(filename)

def test_moldeninput(filepath_molden):
    delete_database('initconds.db')
    initconds = InitialConditionsFactory.from_inputfile(filepath_molden)
    assert initconds.get_condition(50).crd.shape == (10,3)
    assert initconds.get_condition(51).veloc.shape == (10,3)
    assert np.array_equal(initconds.equilibrium.crd, pyrazine_crd) 

def test_freqinput(filepath_freq):
    delete_database('initconds.db')
    initconds = InitialConditionsFactory.from_inputfile(filepath_freq)
    assert np.array_equal(initconds.get_condition(0).crd, [0., 0., 0.])
    assert initconds.get_condition(10).veloc.size == 3
