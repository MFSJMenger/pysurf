import os
import matplotlib.pyplot as plt
import numpy as np

from pysurf.sampling import Sampling
from pysurf.colt import Colt
from pysurf.qctools.converter import Converter, energy_converter

class PlotSpectrum(Colt):
    specfile = os.path.join('spectrum', 'spectrum.db') 
    
    _questions="""
    units = eV :: str :: [au, eV]
    broadening = Lorentzian :: str :: [Gaussian, Lorentzian, None]
    width = 0.1 :: float
    energy_start = 0.0 :: float
    energy_end = 10.0 :: float
    plot_points = 100 :: int
    show_plot = yes :: str :: [yes, no]
    """

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        self.config = config
        sampling = Sampling.from_db(self.specfile)

        nstates = sampling.info['dimensions']['nstates']
        npoints = sampling.nconditions
        data = []
        for  energy, fosc in  zip(sampling._db['energy'], sampling._db['fosc']):
            for idx, en in enumerate(energy[1:]):
                data += [[en - energy[0], fosc[idx+1]]]
        data = np.array(data)

        #change units:
        converter = energy_converter.get_converter(tin='au', tout=config['units'])
        data[:, 0] = converter(data[:, 0])
        if config['broadening'] != None:
            cont_data = self.broadening(data)
        else:
            cont_data = np.sort(data)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title('Spectrum')
        ax.plot(cont_data[:, 0], cont_data[:, 1], color='red')
        plt.show()

    def broadening(self, data):
        res = np.zeros((self.config['plot_points'], 2))
        res[:,0] = np.linspace(self.config['energy_start'], self.config['energy_end'], self.config['plot_points'])
        for  (e0, fosc) in data:
            if self.config['broadening'] == 'Lorentzian':
                func = l_vec(fosc, e0, self.config['width']) 
            if self.config['broadening'] == 'Gaussian':
                func = g_vec(fosc, e0, self.config['width'])
            res[:, 1] += func(res[:, 0])
        res[:, 1] = res[:, 1] / len(data)
        return res
            
def l_vec(fosc, e0, w):
    return np.vectorize(lambda e:fosc * 1 / (1 + ((e - e0)/w*2)**2))

def g_vec(fosc, e0, w):
    return np.vecotrize(lambda e: fosc * np.exp(-np.ln(2)*((e - e0)/w*2)**2))


if __name__ == "__main__":
    PlotSpectrum.from_commandline()
