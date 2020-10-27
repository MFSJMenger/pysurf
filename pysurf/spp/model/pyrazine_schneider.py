import numpy as np
from copy import deepcopy
from ..methodbase import Model
from ...system import Mode



ev2au = 1./27.2113961

class PyrazineSchneider(Model):

    implemented = ["energy", "gradient", "fosc"]
    frequencies = np.array([0.126*ev2au, 0.074*ev2au, 0.118*ev2au])
    masses = np.array(np.ones(3)/frequencies)
    displacements = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    modes = [Mode(freq, dis) for freq, dis in zip(frequencies, displacements)]
    crd = [0.0, 0.0, 0.0]
    nstates = 3

    @classmethod
    def from_config(cls, config):
        return cls()

    def __init__(self):
        """Pyrazine model according to Schneider, Domcke, Koeppel; JCP 92, 1045,
           1990 with 3 modes and 2 excited states.
           The model is in dimensionless normal modes.
        """
        ev2au = 1./27.2113961
        self.E1 = 3.94*ev2au
        self.E2 = 4.84*ev2au
        self.kappa11 = 0.037*ev2au
        self.kappa12 = -0.105*ev2au
        self.kappa21 = -0.254*ev2au
        self.kappa22 = 0.149*ev2au
        self.lam = 0.262*ev2au
        self.w1 = self.frequencies[0]
        self.w2 = self.frequencies[1]
        self.w3 = self.frequencies[2]
        pass

    def get(self, request):
        """the get function returns the adiabatic energies as well as the
           gradient at the given position crd. Additionally the masses
           of the normal modes are returned for the kinetic Hamiltonian.
        """
        crd = request.crd
        T = None
        if 'energy' in request.keys():
            en, T = self.adiab_en(crd)
            request.set('energy', en)
        if 'fosc' in request.keys():
            request.set('fosc', self.adiab_fosc(crd, T))
        if 'gradient' in request.keys():
            grad = self.adiab_grad(crd)
            request.set('gradient', grad)
        return request

    def adiab_en(self, crd):
        """adiab_en returns a one dimensional vector with the adiabatic energies at
           the position crd
        """
        diab_en = self.diab_en(crd)
        (adiab_en, T) = np.linalg.eigh(diab_en)
        adiab_en = np.sort(adiab_en)
        return adiab_en, T

    def adiab_fosc(self, crd,  T=None):
        diab_fosc = [0, 0, 1]
        if T is None:
            diab_en = self.diab_en(crd)
            (en, T) = np.linalg.eigh(diab_en)
        return np.abs(T.dot(diab_fosc))

    def diab_en(self, crd):
        """diab_en returns the full diabatic matrix of the system
        """
        diab_en = np.empty((3, 3), dtype=float)
        diab_en[0, 0] = 0.5*(self.w1*crd[0]**2
                             + self.w2*crd[1]**2
                             + self.w3*crd[2]**2)
        diab_en[0, 1] = 0.0
        diab_en[0, 2] = 0.0
        diab_en[1, 1] = diab_en[0, 0] + self.E1 \
            + self.kappa11*crd[0] + self.kappa12*crd[1]
        diab_en[1, 2] = self.lam*crd[2]
        diab_en[2, 2] = diab_en[0, 0] + self.E2 \
            + self.kappa21*crd[0] + self.kappa22*crd[1]
        diab_en[1, 0] = diab_en[0, 1]
        diab_en[2, 0] = diab_en[0, 2]
        diab_en[2, 1] = diab_en[1, 2]
        return diab_en

    def adiab_grad(self, crd):
        """adiab_grad calculates and returns the numeric adiabatic gradients
        """
        adiab_grad = np.empty((3, 3), dtype=float)
        dq = 0.01
        for i in range(3):
            crd1 = deepcopy(crd)
            crd1[i] = crd1[i] + dq
            en1, _ = self.adiab_en(crd1)
            crd2 = deepcopy(crd)
            crd2[i] = crd2[i] - dq
            en2, _ = self.adiab_en(crd2)
            adiab_grad[:, i] = (en1 - en2)/2.0/dq
        return {i: adiab_grad[i] for i in range(3)}

    def diab_grad(self, crd):
        """diab_grad returns the matrix with the analytical diabatic gradients
        """
        diab_grad = np.empty((3, 3, 3), dtype=float)
        diab_grad[0, 0, 0] = self.w1*crd[0]
        diab_grad[0, 0, 1] = self.w2*crd[1]
        diab_grad[0, 0, 2] = self.w3*crd[2]
        diab_grad[0, 1, 0:3] = 0.0
        diab_grad[0, 2, 0:3] = 0.0
        diab_grad[1, 1, 0] = diab_grad[0, 0, 0] + 2.*self.kappa11*crd[0]
        diab_grad[1, 1, 1] = diab_grad[0, 0, 1] + 2.*self.kappa12*crd[1]
        diab_grad[1, 1, 2] = diab_grad[0, 0, 2]
        diab_grad[1, 2, 0:2] = 0.0
        diab_grad[1, 2, 2] = self.lam
        diab_grad[2, 2, 0] = diab_grad[0, 0, 0] + 2.*self.kappa21*crd[0]
        diab_grad[2, 2, 1] = diab_grad[0, 0, 1] + 2.*self.kappa12*crd[1]
        diab_grad[2, 2, 2] = diab_grad[0, 0, 2]
        diab_grad[1, 0, :] = diab_grad[0, 1, :]
        diab_grad[2, 0, :] = diab_grad[0, 2, :]
        diab_grad[2, 1, :] = diab_grad[1, 2, :]
        return diab_grad


#class PyrMiniSala(PyrMini):
#
#    def __init__(self):
#        """ Pyrazine model according to Sala, Lasorne, Gatti and Guerin;
#            PCCP 2014, 16, 15957 with 3 modes and 2 excited states.
#            The model is in dimensionless normal modes.
#        """
#        ev2au = 1./27.2113961
#        self.E1 = 3.93*ev2au
#        self.E2 = 4.79*ev2au
#        self.kappa11 = -0.038*ev2au
#        self.kappa12 = -0.081*ev2au
#        self.kappa21 = -0.183*ev2au
#        self.kappa22 = 0.128*ev2au
#        self.lam = 0.195*ev2au
#        self.w1 = 0.126*ev2au
#        self.w2 = 0.074*ev2au
#        self.w3 = 0.116*ev2au
#        self.masses = np.array((1.0/self.w1, 1.0/self.w2, 1.0/self.w3))
#        pass


if __name__ == "__main__":
    """small test function printing the energies and gradients at the equilibrium
       point
    """
    crd = np.array([0.0, 0.0, 0.0])
    model = PyrMini()
    print(model.w1, model.w2, model.w3)
    # print(model.get(crd))
