import os
import numpy as np
import matplotlib.pyplot as plt

from pysurf.utils.design import colors
from pysurf.colt import FromCommandline

@FromCommandline("""
infile = population.dat :: file_exists
outfile = population.png :: file
""")
def plot_population(infile, outfile):

    if os.path.isfile(infile):
        with open(infile) as infilev:
            lines = infilev.readlines()
        nstates = len(lines[0].split())-1
        nsteps = len(lines)
        x = np.empty(nsteps, dtype=float)
        y = np.empty((nsteps, nstates), dtype=float)
        
        vfloat = np.vectorize(float)
        for i, line in enumerate(lines):
            y[i,:] = vfloat(line.split()[1::])
            x[i] = float(line.split()[0])

        black = (0,0,0)                                                         
        blue = (63/255, 81/255, 181/255)                                        
        red = (253/255, 86/255, 33/255)
        yellow = (255/255, 233/255, 75/255)
        green = (138/255, 193/255, 73/255)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('steps')
        ax.set_ylabel('population')
        for i in range(nstates):
            ax.plot(x, y[:,i], color=next(colors), label=f'S{i}')
        ax.legend()
        plt.savefig(outfile)
        plt.show()

    else:
        print('Error: file does not exist')
        exit()

if __name__=="__main__":
    plot_population()
