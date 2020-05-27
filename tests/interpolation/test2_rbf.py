import os
import numpy as np
import matplotlib.pyplot as plt

from pysurf.utils import exists_and_isfile
from pysurf.spp import RbfInterpolator
from pysurf.spp.model import PyrazineSchneider, PyrazineSala, HarmonicOscillator
from pysurf.spp import Request
from pysurf.database import PySurfDB
from pysurf.logger import Logger, get_logger

def fill_db(filename, npoints):
    if exists_and_isfile(filename) is True:
        os.remove(filename)


    ho = HarmonicOscillator()
    model = PyrazineSala({'n_states': 3})
    nmodes = len(model.crd)
    nstates = 3
    db = PySurfDB.generate_database(filename, data=['crd', 'energy'], dimensions={'nmodes':nmodes, 'nstates': nstates}, model=True)
    x = (np.random.random((npoints, nmodes))-0.5) * 20
    y = []
    for r in x:
        db.append('crd', r)
        y += [model.get(Request(r, ['energy'], [i for i in range(nstates)]))['energy']]
        db.append('energy', model.get(Request(r, ['energy'], [i for i in range(nstates)]))['energy'])
        db.increase
    return db

#db1 = fill_db('db1.dat', 5000)
db1 = PySurfDB.load_database('db.dat', read_only=True)
#db2 = fill_db('db2.dat', 1000)
db2 = PySurfDB.load_database('prop.db', read_only=True)

nstates=4
nmodes = 9
logger = get_logger('test.log', 'test')
config = {'energy_threshold': 0.02, 'trust_radius_ci': 0.5, 'trust_radius_general': 1.0, 'inverse_distance': False}
rbf = RbfInterpolator(config, db1, ['energy'], nstates, nmodes, logger=logger)

counter = 0
diff_sum = 0
diff = []
for crd, energy in zip(db2['crd'], db2['energy']):
    counter += 1
    res, trust = rbf.get(Request(crd, ['energy'], [0, 1, 2]))
    diff += [res['energy'] - energy]
    diff_sum += (res['energy'] - energy)**2
print('RMS [eV]= ', np.sqrt(diff_sum/counter)*27.2114)
print('Max diff [eV]=', np.amax(diff, axis=0)*27.2114)


