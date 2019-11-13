import numpy as np
from abc import abstractmethod, ABC


class Model(ABC):

    @abstractmethod
    def get(self):
        pass


class TestDB(Model):
    """ class to test the interpolation routine of the database.
        It provides two simple potential energy surfaces and the
        corresponding gradients of shifted harmonic oscillators in an
        arbitrary dimension.
    """
    def __init__(self):
        pass

    def get(self, request):
        """ the get function returns the adiabatic energies as well as
            the gradient at the given position crd. Additionally the
            masses of the normal modes are returned for the kinetic
            Hamiltonian.
        """
        if 'crd' in request.keys():
            crd = request['crd']
            if 'energy' in request.keys():
                request['energy'] = self.en(crd)
            if 'gradient' in request.keys():
                request['gradient'] = self.grad(crd)

        # If an empty dictionary is sent, all possible properties are
        # sent back
        else:
            request['energy'] = None
            request['gradient'] = None
        return request
            

    def en(self, crd):
        """ Method to calculate the energies of the PE surfaces at the
            given point.
        """
        en = [0.0, 0.0]
        for co in crd.flatten():
            en[0] += co**2
            en[1] += (co-1)**2
        return en

    def grad(self, crd):
        """ Method to calculate the gradient for the PE surfaces
            at the given point.
        """
        grad = np.zeros((2, *crd.shape), dtype=float)
        grad[0, :, :] = 2.*crd
        grad[1, :, :] = 2.*(crd-1)
        return grad


if __name__ == "__main__":
    """ small test function printing the energies and gradients at the
        equilibrium point
    """
    crd = np.array([[0.0, 0.0, 0.0]])
    model = TestDB()
    print(model.get(crd))
