import os
import numpy as np


listdir = os.listdir('./')
vfloat = np.vectorize(float)

atoms = ['C','N','C','C','H','H','N','C','H','H']
bohr2angstrom = 0.529177208

for obj in listdir:
    if obj.startswith('traj.'):
        if os.path.isfile(obj + '/crd.txt'):
            with open(obj + '/crd.txt') as infile:
                inp = infile.readlines()
            with open(obj + '/crd.xyz', 'w') as outfile:
                for i in range(1,len(inp)):
                    outfile.write('10\n')
                    outfile.write(str(i)+'\n')
                    crds = vfloat(inp[i].split()).reshape(10,3)
                    for j in range(len(crds)):
                        outfile.write("{0}  {1:12.8f}  {2:12.8f}  {3:12.8f}\n".format(atoms[j], crds[j,0]*bohr2angstrom,
                        crds[j,1]*bohr2angstrom, crds[j,2]*bohr2angstrom))
