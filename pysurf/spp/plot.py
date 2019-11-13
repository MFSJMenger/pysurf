import matplotlib.pyplot as plt
import os
import sys
import numpy as np

from spp import SurfacePointProvider

if os.path.abspath('../') not in sys.path:
    sys.path.insert(0, os.path.abspath('../'))
from utils.strutils import split_str


class PlotPES(SurfacePointProvider):
    """ Class to plot the PES of the SPP
    """
    def __init__(self, inputfile):
        """ Here, the same logger is used as in SPP.
            Here, the same input parser is used  as in SPP
        """
        super(PlotPES, self).__init__(inputfile)

        #check whether there is a PLOT section in the inputfile
        if 'PLOT' not in self.config.keys():
            self.logger.error('No plot section in inputfile!')
            exit()

        #check whether the inputfile contains which coordinate is used
        #for the plotting
        if 'crd' not in self.config['PLOT'].keys():
            self.logger.error('No coordinate given in plot section '
                              + 'of the inputfile')

        if len(split_str(self.config['PLOT']['crd'])) == 1:
            try:
                self.plot_crd = int(self.config['PLOT']['crd'])
            except ValueError:
                self.logger.error('crd in PLOT section is not an integer!')
                exit()
            #start plotting of 1D PES
            self.plot1DPES()
        else:
            self.logger.info('plots in higher dimensions are not '
                             + 'yet implemented')

    def read_range(self):
        """ method to read the range input from the PLOT section
            of the input file. Error handling is included.
        """
        if 'range' in self.config['PLOT'].keys():
            try:
                ran = split_str(self.config['PLOT']['range'])
                xrange = [float(st) for st in ran]
            except (KeyError,ValueError):
                self.logger.error('Range input not correct')
                exit()
            if len(xrange) != 2:
                self.logger.error('Range input not correct')
                exit()
        else:
            xrange = (-5, 5)
        return xrange

    def read_npoints(self):
        """ method to read the plotting grid from the inputfile.
            error handling is included.
        """
        if 'number of points to plot' in self.config['PLOT'].keys():
            try:
                npoints = int(self.config['PLOT']['number of points to plot'])
            except ValueError:
                self.logger.warning('Number of points to plot argument is not '
                                    + 'an integer! Using default!')
                npoints = 100
        else:
            npoints = 100
        return npoints

    def read_nstates(self):
        """ method to get the states for which the PESs have to be 
            plotted. Error handling is included.
        """
        if 'number of states to plot' in self.config['PLOT'].keys():
            try:
                nstates = int(self.config['PLOT']['number of states to plot'])
            except ValueError:
                self.logger.warning('Number of states to plot argument is not '
                                    + 'an integer! Using default!')
                nstates = 1
        else:
            nstates = 1
        return nstates

    def read_ndofs(self):
        """ method to read the number of degrees of the system from the input
            file. Error handling is included.
        """
        if 'number of degrees of freedom' not in self.config['PLOT'].keys():
            self.logging.error('Please provide number of degrees of freedom '
                               + ' in the plot section!')
            exit()
        else:
            try:
                ndofs = int(self.config['PLOT']['number of degrees '
                                                + 'of freedom'])
            except ValueError:
                self.logger.error('Number of degrees of freedom is not'
                                  + ' an integer!')
                exit()
        return ndofs

    def plot1DPES(self):
        """ method to create a 1D plot of the PES according to 
            the input parameters from the input file
        """
        xrange = self.read_range()
        nstates = self.read_nstates()
        npoints = self.read_npoints()
        ndofs = self.read_ndofs()
        xval = np.linspace(*xrange, npoints)
        yval = np.empty((npoints, nstates), dtype=float)

        for i in range(0, len(xval)):
            crd = np.zeros(ndofs)
            crd[self.plot_crd-1] = xval[i]
            yval[i, :] = self.get(crd)['energy']

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel('Coord ' + str(self.plot_crd))
        ax.set_ylabel('energy')

        ax.plot(xval, yval)
        plt.show()


if __name__ == '__main__':
    myplot = PlotPES('test.inp')
