import numpy as np
from copy import deepcopy
from ..methodbase import ModelBase
from pysurf.system import Mode


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
        self.a = 1.0
        self.x0 = 0.0

    def energy(self, x):
        return 0.5*self.a*(x - self.x0)**2

    def gradient(self, x):
        return self.a*(x - self.x0)

    def get(self, request):
        """the get function returns the adiabatic energies as well as the
           gradient at the given position crd. Additionally the masses
           of the normal modes are returned for the kinetic Hamiltonian.
        """
        crd = request['crd']
        en = self.energy(crd)
        grad = self.gradient(crd)
        return {'crd': crd, 'energy': en, 'gradient': {0: grad}, 'masses': self.masses}
