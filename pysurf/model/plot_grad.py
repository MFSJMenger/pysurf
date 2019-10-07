import numpy as np
import matplotlib.pyplot as plt

from pyrmini import PyrMini

x = np.linspace(-5,5,100)

y_grad = np.empty((100,3),dtype=float)
y_en = np.empty((100,3), dtype=float)

pm=PyrMini()
print(pm.w1)
for i in range(100):
    res = pm.get(np.array([0,0,x[i]]))
    y_grad[i] = res['gradient'][:,2]
    y_en[i] = res['energy']
    #print(pm.get(np.array([x[i],0,0]))['gradient'])


fig = plt.figure()
ax = fig.add_subplot(1, 2, 1)
ax.set_xlabel('Coord ')
ax.set_ylabel('energy')
ax.plot(x, y_grad)
ax = fig.add_subplot(1, 2, 2)
ax.set_xlabel('Coord ')
ax.set_ylabel('energy')
ax.plot(x, y_en)
plt.show()

