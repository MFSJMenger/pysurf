import numpy as np
from collections import namedtuple

from pysurf.spp import Model
from pysurf.system import Mode


class HarmonicOscillator1D(Model):
    """ Model for a 1D harmonic oscillator with 1 or 2 potential energy surfaces """

    _questions = """
    e0 = 0.0 :: float
    w0 = 1.0 :: float
    x0 = 0.0 :: float
    # Number of potential energy surfaces
    npes = 1 :: str ::

    [npes(1)]

    [npes(2)]
    e1 = 1.0 :: float
    w1 = 1.0 :: float
    x1 = 1.0 :: float 
    """

    implemented = ["energy", "gradient"]
    masses = [1.0]
    crd = [0.0]
    frequencies = [1.0]
    displacements = [[1.0]]
    modes = [Mode(freq, dis) for freq, dis in zip(frequencies, displacements)]

    @classmethod
    def from_config(cls, config):
        e0 = config['e0']
        w0 = config['w0']
        x0 = config['x0']
        npes = int(config['npes'].value)
        #
        config_npes = config['npes']
        return cls(e0, w0, x0, npes, config_npes)

    def __init__(self, e0, w0, x0, npes, config_npes):
        """ init function for the 1D Harmonic Oscillator. Take care that the mass is set to 1, i.e.
            fs are not an appropriate timescale for this model, but atomic units!

            Parameters:
            -----------
                e0: float
                    Float number for the energy offset of the lowest PE surface

                w0: float
                    Float number for the frequency of the lowest PE surface

                x0: float
                    Float number for the shift of the lowest PE surface

                npes: int
                    1 or 2, depending whether 1 or 2 PE surfaces should be in the model

                config_npes:
                    config object (similar to dictionary) containing the information of the 
                    second PE surface. If npes == 1, it is not used.
        """
                    
        self.frequencies = [w0]
        self.crd = [x0]
        self.npes = int(npes)
        self.w = [w0]
        self.x = [x0]
        self.e = [e0]
        #
        if self.npes == 2:
            self.w += [config_npes['w1']]
            self.x += [config_npes['x1']]
            self.e += [config_npes['e1']]


    def _energy(self, x):
        """ returns the energies of the model at position x """
        energy = []
        for i in range(self.npes):
            energy += [0.5*self.w[i]*(x - self.x[i])**2 + self.e[i]]
        energy = np.array(energy).flatten()
        return energy

    def _gradient(self, x):
        """ returns the gradients of the model at position x """
        gradient = {}
        for i in range(self.npes):
            gradient[i] = np.array(self.w[i]*(x - self.x[i]))
        return gradient
         

    def get(self, request):
        """ the get function returns the adiabatic energies as well as the
            gradient at the given position request.crd.
           
            Parameters:
            -----------
                request:
                    request object containing the information, which properties are asked for


            Returns:
            --------
                request:
                    request object with the asked information
        """
        crd = request.crd
        for prop in request:
            if prop == 'energy':
                request.set('energy', self._energy(crd))
            if prop == 'gradient':
                request.set('gradient', self._gradient(crd))
        return request


if __name__=="__main__":
    HO = HarmonicOscillator1D(E0=0.0, w0=1.0, x0=0.0, npes=1, config_npes={})
    fake_request = namedtuple('request', 'crd energy gradient')
    fake_request.crd = np.array([1.0])
    print(HO._energy(fake_request.crd))
