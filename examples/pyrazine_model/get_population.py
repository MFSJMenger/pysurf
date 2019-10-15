import os
import numpy as np

lsdir = os.listdir('./')

counter = 0
population = np.zeros((10001,3), dtype=int)
for obj in lsdir:
    if obj.startswith('traj.'):
        counter += 1
        with open(obj + '/energy.txt') as infile:
            inp = infile.readlines()
        for i in range(1, len(inp)):
            state = int(inp[i].split()[0])
            population[i,state] += 1

population = population/float(counter)

np.savetxt('pop.dat',population)


