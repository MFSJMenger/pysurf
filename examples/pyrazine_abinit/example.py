import numpy as np
import os

from pysurf.wigner import InitialCondition
from pysurf.wigner import get_initial_condition
from pysurf.wigner import Molecule
from pysurf.wigner import Mode
from pysurf.wigner import WignerSampling

from pysurf.model.pyrmini import PyrMini
import numpy.random as random

from model import landau_zener_surfacehopping


def get_init(init=0, seed=16661):
    random.seed(seed)
    
    wigner = WignerSampling.from_molden('/data/ehrmaier/pysurf/examples/pyrazine_abinit/preparation/molden.in')
    
    return wigner.create_initial_conditions('/data/ehrmaier/pysurf/examples/pyrazine_abinit/preparation/molden.in',init)


init = get_init(1)
print(init)
# create initial condition
#for i in range(100):
#    os.mkdir('traj.'+str(i))
#    os.chdir('traj.'+str(i))
#    init = get_init(i)
#    os.system('cp ../test.inp ./')
#    landau_zener_surfacehopping(init, 2, 10000, 16661, 'test.inp', 10)
#    os.chdir('../')

