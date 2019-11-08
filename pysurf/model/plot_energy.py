import numpy as np
import matplotlib.pyplot as plt

from pyrmini import PyrMini

x = np.linspace(-10,10,100)

y_grad = np.empty((100,3),dtype=float)
y1_en = np.empty((100,3), dtype=float)
y2_en = np.empty((100,3), dtype=float)
y3_en = np.empty((100,3), dtype=float)

pm=PyrMini()
print(pm.w1)
for i in range(100):
    res = pm.get({'coord': np.array([x[i],0,0]), 'gradient': None, 'energy': None})
    y1_en[i] = res['energy']
    #print(pm.get(np.array([x[i],0,0]))['gradient'])

for i in range(100):
    res = pm.get({'coord': np.array([0,x[i],0]), 'gradient': None, 'energy': None})
    y2_en[i] = res['energy']
    #print(pm.get(np.array([x[i],0,0]))['gradient'])

for i in range(100):
    res = pm.get({'coord': np.array([0,0,x[i]]), 'gradient': None, 'energy': None})
    y3_en[i] = res['energy']

y1_en = np.array(y1_en)*27.2114
y2_en = np.array(y2_en)*27.2114
y3_en = np.array(y3_en)*27.2114


black = (0,0,0)
blue = (63/255, 81/255, 181/255)
red = (253/255, 86/255, 33/255)
fig = plt.figure()
ax = fig.add_subplot(1, 3, 1)
ax.set_xlabel(r'$Q_1$')
ax.set_ylabel('energy [eV]')
for i, color in enumerate((black, blue, red)):
    ax.plot(x, y1_en[:,i], color=color)
ax = fig.add_subplot(1, 3, 2)
ax.set_xlabel(r'$Q_{6a}$')
ax.set_ylabel('energy [eV]')
for i, color in enumerate((black, blue, red)):
    ax.plot(x, y2_en[:,i], color=color)
ax = fig.add_subplot(1, 3, 3)
ax.set_xlabel(r'$Q_{10a}$')
ax.set_ylabel('energy [eV]')
for i, color in enumerate((black, blue, red)):
    ax.plot(x, y3_en[:,i], color=color)
plt.show()

