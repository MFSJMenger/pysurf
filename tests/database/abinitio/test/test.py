from scipy.interpolate import LinearNDInterpolator
import numpy as np

x = np.linspace(0,1,100)

points = np.empty((10000,2), dtype=float)
z = np.empty((10000,3), dtype=float)

for i in range(100):
    for j in range(100):
        z[(j*i) + j,0] = (i/100)**2+(j/100.)**2
        z[(j*i) + j,1] = 1
        z[(j*i) + j,2] = 2 
        points[(j*i) + j,:] = [i/100.,j/100.]

print(points)
print(z)
interpolator = LinearNDInterpolator(points, z)

print(interpolator((0.24,0.222)))

    
