import numpy as np
import os
import shutil
from copy import deepcopy

from pysurf.wigner import InitialCondition
from pysurf.wigner import get_initial_condition
from pysurf.wigner import Molecule
from pysurf.wigner import Mode
from pysurf.wigner import WignerSampling
from pysurf.spp.spp import SurfacePointProvider

#from pysurf.model.pyrmini import PyrMini
import numpy.random as random

from pysurf.sh.model import landau_zener_surfacehopping


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
ntrajs = 20
# create initial condition
conditions = get_init(ntrajs)
#print(conditions.get_condition(0))

#spp = SurfacePointProvider('test_db.inp')
#refgeo = conditions.equilibrium.crd
#print('Johannes refgeo: ',refgeo)
#coords = deepcopy(refgeo)
#coords[1,0] += 0.3
#coords[2,0] +=0.3
#res = spp.get({'coord': coords, 'energy': None, 'gradient': None})
#coords[3,0] += 0.3
#res = spp.get({'coord': coords, 'energy': None, 'gradient': None})
#coords = deepcopy(refgeo)
#coords[1,0] += 0.1
#res = spp.get({'coord': coords, 'energy': None, 'gradient': None})
#for i in range(len(refgeo)):
#    for j in range(len(refgeo[i])):
#        coords = deepcopy(refgeo)
#        coords[i,j] += 0.1
#        res = spp.get({'coord': coords, 'energy': None, 'gradient': None})
#        print(res)

for i in range(ntrajs):
    os.mkdir('traj.'+str(i))
    os.chdir('traj.'+str(i))
    os.system('cp ../test_db.inp ../ref_geo.xyz ../template.dat ../qchem.out ../run_qchem ../db.dat ./')
    landau_zener_surfacehopping(conditions.get_condition(i), 2, 200, 16661, 'test_db.inp', 20)
    os.system('cp db.dat ../')
    os.chdir('../')

