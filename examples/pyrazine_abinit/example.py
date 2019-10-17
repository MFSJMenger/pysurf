import numpy as np
import os
import shutil

from pysurf.wigner import InitialCondition
from pysurf.wigner import get_initial_condition
from pysurf.wigner import Molecule
from pysurf.wigner import Mode
from pysurf.wigner import WignerSampling

#from pysurf.model.pyrmini import PyrMini
import numpy.random as random

from model import landau_zener_surfacehopping


print('Johannes hallo!')

def get_init(init=0, seed=16666):
    random.seed(seed)
    # set up wigner sampling from molden file    
    sampling = WignerSampling.from_molden('/data/ehrmaier/pysurf/examples/pyrazine_abinit/preparation/molden.in')
    
    # create an initial conditions file
    conditions = sampling.create_initial_conditions('init_conds.db',init)
    return conditions



if os.path.isfile('init_conds.db'):
    os.remove('init_conds.db')
if os.path.isdir('traj.0'):
    shutil.rmtree('traj.0')
if os.path.isdir('traj.1'):
    shutil.rmtree('traj.1')
# create initial condition
conditions = get_init(2)
print(conditions.get_condition(0))
for i in range(2):
    os.mkdir('traj.'+str(i))
    os.chdir('traj.'+str(i))
    os.system('cp ../test_abinit.inp ../ref_geo.xyz ../template.dat ../qchem.out ../run_qchem ./')
    landau_zener_surfacehopping(conditions.get_condition(i), 2, 100, 16661, 'test_abinit.inp', 10)
    os.chdir('../')

