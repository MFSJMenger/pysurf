import numpy as np
import random
import matplotlib.pyplot as plt

from pysurf.spp.dbinter.dbinter import ShepardInterpolator
from pysurf.spp.dbinter.dbinter import MyRBF

start = -1.2
stop = 1.2

x = np.array([(random.random()-0.5)*2. for _ in range(10)])
print(x)
y = x**2
print(y)


interpolator = ShepardInterpolator(x, y)
rbfinterpolator = MyRBF(x, y)
print(rbfinterpolator)

x_inter = np.linspace(start,stop,50)
y_inter = np.array([interpolator(i) for i in x_inter]) 

y_rbf = np.array([rbfinterpolator([i]) for i in x_inter])

x_real = np.linspace(start, stop, 50)
y_real = x_real**2

black = (0,0,0)                                                         
blue = (63/255, 81/255, 181/255)                                        
red = (253/255, 86/255, 33/255)   

fig = plt.figure()
ax = fig.add_subplot(1,2,1)
ax.plot(x_real, y_real, linestyle= 'dashed',color='black')
ax.plot(x_inter, y_inter, color=blue)
ax.plot(x, y, 'o', color=red)

ax = fig.add_subplot(1,2,2)
ax.plot(x_real, y_real, linestyle= 'dashed',color='black')
ax.plot(x_inter, y_rbf, color=blue)
ax.plot(x, y, 'o', color=red)






plt.show()



