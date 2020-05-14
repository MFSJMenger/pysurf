import os
import numpy as np
import matplotlib.pyplot as plt

from pysurf.utils import exists_and_isfile
from pysurf.spp import RbfInterpolator
from pysurf.spp.model import PyrMini, HarmonicOscillator
from pysurf.spp import Request
from pysurf.database import PySurfDB


def fill_db(filename, npoints):
    if exists_and_isfile(filename) is True:
        os.remove(filename)
    db = PySurfDB.generate_database(filename, data=['crd', 'energy'], dimensions={'nmodes':3, 'nstates': 2}, model=True)
    x = (np.random.random(npoints)-0.5) * 10
    x2 = (np.random.random(npoints)-0.5)*10
    x3 = (np.random.random(npoints)-0.5)*10

    ho = HarmonicOscillator()
    pyr = PyrMini()
    y = []
    for r1, r2, r3 in zip(x, x2, x3):
        db.append('crd', [r1, r2, r3])
        y += [pyr.get(Request([r1, r2, r3], ['energy'], [1, 2]))['energy'][1:]]
        db.append('energy', pyr.get(Request([r1, r2, r3], ['energy'], [1, 2]))['energy'][1:])
        db.increase
    return db

db1 = fill_db('db1.dat', 1000)
db2 = fill_db('db2.dat', 1000)

rbf = RbfInterpolator(db1, ['energy'])

counter = 0
diff = 0
for crd, energy in zip(db2['crd'], db2['energy']):
    counter += 1
    res, trust = rbf.get(Request(crd, ['energy'], [1, 2]))
    diff += (res['energy'] - energy)**2
print('RMS [eV]= ', np.sqrt(diff/counter)*27.2114)


