import numpy as np
from copy import deepcopy
from abc import abstractmethod, ABC 

class Model(ABC):

    @abstractmethod
    def get(self):
        pass



class PyrMini(Model):

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
        self.kappa21 = -0.0254*ev2au
        self.kappa22 = 0.149*ev2au
        self.lam = 0.262*ev2au
        self.w1 = 0.126*ev2au
        self.w2 = 0.074*ev2au
        self.w3 = 0.118*ev2au
        self.mass = np.array((1.0/self.w1, 1.0/self.w2, 1.0/self.w3))
        pass

    def get(self, request):
        """the get function returns the adiabatic energies as well as the
           gradient at the given position coord. Additionally the masses
           of the normal modes are returned for the kinetic Hamiltonian.
        """
        if 'coord' in request.keys():
            coord = request['coord']
            if 'energy' in request.keys():
                request['energy'] = self.adiab_en(coord)
            if 'gradient' in request.keys():
                request['gradient'] = self.adiab_grad(coord)
            request['mass'] = self.mass
        else:
            request['energy'] = None
            request['gradient'] = None
            request['mass'] = None
        return request

    def adiab_en(self, coord):
        """adiab_en returns a one dimensional vector with the adiabatic energies at
           the position coord
        """
        diab_en = self.diab_en(coord)
        adiab_en = np.linalg.eig(diab_en)[0]
        adiab_en = np.sort(adiab_en)
        return adiab_en

    def diab_en(self, coord):
        """diab_en returns the full diabatic matrix of the system
        """
        diab_en = np.empty((3, 3), dtype=float)
        diab_en[0, 0] = 0.5*(self.w1*coord[0]**2
                             + self.w2*coord[1]**2
                             + self.w3*coord[2]**2)
        diab_en[0, 1] = 0.0
        diab_en[0, 2] = 0.0
        diab_en[1, 1] = diab_en[0, 0] + self.E1 \
            + self.kappa11*coord[0] + self.kappa12*coord[1]
        diab_en[1, 2] = self.lam*coord[2]
        diab_en[2, 2] = diab_en[0, 0] + self.E2 \
            + self.kappa21*coord[0] + self.kappa22*coord[1]
        diab_en[1, 0] = diab_en[0, 1]
        diab_en[2, 0] = diab_en[0, 2]
        diab_en[2, 1] = diab_en[1, 2]
        return diab_en

    def adiab_grad(self, coord):
        """adiab_grad calculates and returns the numeric adiabatic gradients
        """
        adiab_grad = np.empty((3, 3), dtype=float)
        dq = 0.01
        for i in range(3):
            coord1 = deepcopy(coord)
            coord1[i] = coord1[i] + dq
            en1 = self.adiab_en(coord1)
            coord2 = deepcopy(coord)
            coord2[i] = coord2[i] - dq
            en2 = self.adiab_en(coord2)
            adiab_grad[:, i] = (en1 - en2)/2.0/dq
        return adiab_grad

    def diab_grad(self, coord):
        """diab_grad returns the matrix with the analytical diabatic gradients
        """
        diab_grad = np.empty((3, 3, 3), dtype=float)
        diab_grad[0, 0, 0] = self.w1*coord[0]
        diab_grad[0, 0, 1] = self.w2*coord[1]
        diab_grad[0, 0, 2] = self.w3*coord[2]
        diab_grad[0, 1, 0:3] = 0.0
        diab_grad[0, 2, 0:3] = 0.0
        diab_grad[1, 1, 0] = diab_grad[0, 0, 0] + 2.*self.kappa11*coord[0]
        diab_grad[1, 1, 1] = diab_grad[0, 0, 1] + 2.*self.kappa12*coord[1]
        diab_grad[1, 1, 2] = diab_grad[0, 0, 2]
        diab_grad[1, 2, 0:2] = 0.0
        diab_grad[1, 2, 2] = self.lam
        diab_grad[2, 2, 0] = diab_grad[0, 0, 0] + 2.*self.kappa21*coord[0]
        diab_grad[2, 2, 1] = diab_grad[0, 0, 1] + 2.*self.kappa12*coord[1]
        diab_grad[2, 2, 2] = diab_grad[0, 0, 2]
        diab_grad[1, 0, :] = diab_grad[0, 1, :]
        diab_grad[2, 0, :] = diab_grad[0, 2, :]
        diab_grad[2, 1, :] = diab_grad[1, 2, :]
        return diab_grad


if __name__ == "__main__":
    """small test function printing the energies and gradients at the equilibrium
       point
    """
    coord = np.array([0.0, 0.0, 0.0])
    model = PyrMini()
    print(model.get(coord))
