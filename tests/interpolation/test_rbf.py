import os
import numpy as np
import matplotlib.pyplot as plt

from pysurf.spp import RbfInterpolator
from pysurf.spp.model import PyrMini, HarmonicOscillator
from pysurf.spp import Request
from pysurf.database import PySurfDB


def fill_db():
    os.remove('db.dat')
    db = PySurfDB.generate_database('db.dat', data=['crd', 'energy'], dimensions={'nmodes':3, 'nstates': 2}, model=True)
    x = (np.random.random(5000)-0.5) * 10
    x2 = (np.random.random(5000)-0.5)*10
    x3 = (np.random.random(5000)-0.5)*10

    ho = HarmonicOscillator()
    pyr = PyrMini()
    y = []
    for r1, r2, r3 in zip(x, x2, x3):
        db.append('crd', [r1, r2, r3])
        y += [pyr.get(Request([r1, r2, r3], ['energy'], [1, 2]))['energy'][1:]]
        db.append('energy', pyr.get(Request([r1, r2, r3], ['energy'], [1, 2]))['energy'][1:])
        db.increase

fill_db()
pyr = PyrMini()
x_lin = np.linspace(-5, 5, 10)
y_lin = [pyr.get(Request([0, 0, r], ['energy'], [1, 2]))['energy'][1:] for r in x_lin]

db = PySurfDB.load_database('db.dat', data=['crd', 'energy'], dimensions={'nmodes':3, 'nstates': 2}, model=True)
rbf = RbfInterpolator(db, ['energy'])


x = np.linspace(-10, 10, 100)
y=[]
for r in x:
    request = Request([0, 0, r], ['energy'], [1, 2])
    res, trust = rbf.get(request)
    y += [res['energy']]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(x, np.array(y)[:, 0])
ax.plot(x, np.array(y)[:, 1])
ax.plot(x_lin, y_lin, marker='o', color='red', linestyle='None') 
plt.show()
