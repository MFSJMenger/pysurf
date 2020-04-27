import os
import matplotlib.pyplot as plt
import numpy as np

from pysurf.sampling import Sampling
from pysurf.colt import Colt

class PlotSpectrum(Colt):
    specfile = os.path.join('spectrum', 'spectrum.db') 
    
    _questions="""
    units = au :: str :: [au]
    broadening = Lorentzian :: str :: [Gaussian, Lorentzian]
    width = 0.001 :: float
    energy_start = 0.0 :: float
    energy_end = 0.45 :: float
    plot_points = 100 :: int
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
        if config['broadening'] == 'Lorentzian':
            cont_data = self.lorentzian(data, config['width'])
               
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title('Spectrum')
        ax.stem(data[:,0], data[:,1])
#        ax.plot(cont_data[:,0], cont_data[:,1])
        plt.show()

    def lorentzian(self, data, width):
        res = np.zeros((self.config['plot_points'], 2))
        res[:,0] = np.linspace(self.config['energy_start'], self.config['energy_end'], self.config['plot_points'])
        for  (e0, fosc) in data:
            func = l_vec(fosc, e0, self.config['width']) 
            res[:, 1] += func(res[:, 0])
        res[:, 1] = res[:, 1] / len(data)
        return res
            
def l_vec(fosc, e0, w):
    return np.vectorize(lambda e:fosc * 1 / (1 + ((e - e0)/w*2))**2)


if __name__ == "__main__":
    PlotSpectrum.from_commandline()
