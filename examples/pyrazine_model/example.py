import numpy as np
import os

from pysurf.wigner import InitialCondition
from pysurf.wigner import get_initial_condition
from pysurf.wigner import Molecule
from pysurf.wigner import Mode

from pysurf.model.pyrmini import PyrMini
import numpy.random as random

from model import landau_zener_surfacehopping


def get_init(init=0, seed=16661):
    random.seed(seed)

    obj = PyrMini()
    mol = Molecule(np.array([1, 1, 1]), np.zeros(3, dtype=np.double), np.ones(3, dtype=np.double))
    #
    displacement = np.ones(3, dtype=np.double)
    displacement /= np.linalg.norm(displacement)
    # modes
    mode1 = Mode(obj.w1, displacement)
    mode2 = Mode(obj.w2, displacement)
    mode3 = Mode(obj.w3, displacement)

    data = obj.get({'coord': np.array([1, 1, 1]), 'mass': None, 'energy': None, 'gradient': None})
    mol.masses = data['mass']

    crd = np.array([0., 0., 0.], dtype=np.double)
    veloc = np.array([0., 0., 0.], dtype=np.double)


    for _ in range(init):
        # sample each mode seperately
        for i, mode in enumerate((mode1, mode2, mode3)):
            _, cond = get_initial_condition(mol, [mode])
            crd[i] = cond.crd[i]
            veloc[i] = cond.veloc[0][i]
    return InitialCondition(crd, veloc)
# create initial condition
for i in range(100):
    os.mkdir('traj.'+str(i))
    os.chdir('traj.'+str(i))
    init = get_init(i)
    os.system('cp ../test.inp ./')
    landau_zener_surfacehopping(init, 2, 10000, 16661, 'test.inp', 10)
    os.chdir('../')

