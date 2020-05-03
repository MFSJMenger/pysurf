import numpy as np
from copy import deepcopy
from ..methodbase import ModelBase
from pysurf.molecule import Mode


class HarmonicOscillator(ModelBase):

    implemented = ["energy", "gradient"]
    masses = [1.0]
    crd = [0.0]
    frequencies = [1.0]
    displacements = [[1.0]]
    modes = [Mode(freq, dis) for freq, dis in zip(frequencies, displacements)]

    @classmethod
    def from_config(cls, config):
        return cls()

    def __init__(self):
        self.a = 5.0
        self.x0 = np.array([0.0, 5.0, -2.0])

    def energy(self, x):
        return self.a*(x - self.x0)**2

    def gradient(self, x):
        print(x)
        print(x - self.x0)
        return -2*self.a*(x - self.x0)

    def get(self, request):
        """the get function returns the adiabatic energies as well as the
           gradient at the given position crd. Additionally the masses
           of the normal modes are returned for the kinetic Hamiltonian.
        """
        crd = request['crd']
        val = 1
        N = 3
        en = self.energy(crd)
        grad = self.gradient(crd)
        en = [e + val*i for i, e in enumerate(en)]
        grad = self.gradient(crd)
        grad = [grad for _ in range(N)]
        mass = [1, 1, 1]
        return {'energy': [en], 'gradient': {0: grad}, 'mass': mass}



if __name__ == "__main__":
    """small test function printing the energies and gradients at the equilibrium
       point
    """
    crd = np.array([0.0, 0.0, 0.0])
    model = HarmonicOscillator()
    print(model.w1, model.w2, model.w3)
    # print(model.get(crd))
