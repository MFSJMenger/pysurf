import numpy as np
#
from colt import Colt
from pysurf.system import ATOMID_TO_NAME
from pysurf.logger import Logger, get_logger
#
from pysurf.spp import AbinitioBase


#has to be adapted according to Max's new implementation
from pyscf import gto, dft, tddft, grad

class DFT(Colt):
    """ class which executes the DFT and TDDFT calculations using the PySCF package """

    _user_input = """
    functional = :: str :: ['pbe0']
    """

    implemented = ['energy', 'gradient']


    @classmethod
    def from_config(cls, config, mol, nstates):
        """ """
        functional = config['functional']
        return cls(functional, mol, nstates)
        

    def __init__(self, functional, mol, nstates):
        """
            Parameters:
            -----------
                functional: str
                    string of the functional

                mol:
                    mol object of PySCF

                nstates: int
                    number of states

        """
        self.mol = mol
        self.nstates = nstates
        
        mydft = dft.RKS(mol).x2c().set(xc=functional)
        
        self.dft_scanner = mydft.as_scanner()
        self.dft_grad = mydft.nuc_grad_method().as_scanner()

        if self.nstates > 1:
            # Switch to xcfun because 3rd order GGA functional derivative is not
            # available in libxc
            mydft._numint.libxc = dft.xcfun
            mytddft = tddft.TDDFT(mydft)
            self.tddft_scanner = mytddft.as_scanner()
            self.tddft_scanner.nstates = self.nstates - 1
            
            # PySCF-1.6.1 and newer supports the .Gradients method to create a grad
            # object after grad module was imported. It is equivalent to call the
            # .nuc_grad_method method.
            self.tddft_grad = mytddft.Gradients().as_scanner()


    def do_energy(self, request, mol):
        """ function to calculate the energies 
            Parameters:
            -----------
                request:
                    request object
                mol: 
                    mol object with the correct coordinates
            Return:
            -------
                request where the energies are filled in
        """
        if self.nstates == 1:
            en = [self.dft_scanner(mol)]
        else:
            en = self.tddft_scanner(mol)

        request.set('energy', en)
   
    def do_gradient(self, request, mol):
        """ function to calculate the gradients 
            Parameters:
            -----------
                request:
                    request object
                mol: 
                    mol object with the correct coordinates
            Return:
            -------
                request where the gradients are filled in
        """
        grad = {}
        for state in request.states:
            if state == 0:
                e_gs, grad_gs = self.dft_grad(mol)
                grad[0] = grad_gs
            else:
                e, grad_tddft = self.tddft_grad(mol, state=state)
                grad[state] = grad_tddft
        request.set('gradient', grad)

class PySCF(AbinitioBase):
    """ Interface for the PySCF code, which is free available 

        The communication with the SPP is the get function with the request object.
        The actual calculations are performed in separate classes. The classes need to have
        a function for each property with the name do_prop, e.g. do_energy for energy calculations.
        All properties which are implemented have to be stored in the implemented property of the 
        corresponding classes
    """

    _user_input = """
    basis = 631g*
    # Calculation Method
    method = DFT/TDDFT :: str :: [DFT/TDDFT]
    """

    # implemented has to be overwritten by the individual classes for the methods
    implemented = []

    # dictionary containing the keywords for the method and the corresponding classes
    methods = {'DFT/TDDFT': DFT}

    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                             for name, method in cls.methods.items()})

    @classmethod
    def from_config(cls, config, atomids, nstates, nghost):
        if nghost != 0:
            raise ValueError("Number of ghost states has to be 0 for the PySCF interface")
        method = config['method'].value
        basis = config['basis']
        config_method = config['method']
        return cls(basis, method, atomids, nstates, config_method)


    def __init__(self, basis, method, atomids, nstates, config_method):
        """
            Parameters:
            -----------
                basis: str
                    containing the string of the basis set, directly passed to pyscf

                method: str
                    string of the method, which has to be in the methods dictionary

                atomids: list
                    list with the atomids

                nstates: int
                    number of states, including the ground-state, i.e. 1 is only the ground-state

                config_method:
                    Colt config for the individual method
        """
        self.mol = self._generate_pyscf_mol(basis, atomids)
        self.nstates = nstates
        self.atomids = atomids
        self.basis = basis
        # initializing the class for the corresponding method
        self.calculator = self.methods[method].from_config(config_method, self.mol, nstates)
        # update the implemented property
        self.implemented = self.calculator.implemented


    
    def get(self, request):
        """ main interface with the SPP 

            Paramters:
            ----------
                request:
                    request instance of the Request class containing. It contains the coordinates
                    and the desired properties that should be calculated

            Return:
            -------
                request:
                    request instance of the Request class where the desired information has been
                    filled in
        """
        # update coordinates
        self.mol = self._generate_pyscf_mol(self.basis, self.atomids, request.crd)
        for prop in request:
            func = getattr(self.calculator, 'do_' + prop)
            func(request, self.mol)
        #
        return request

    @staticmethod
    def _generate_pyscf_mol(basis, atomids, crds=None):
        """ helper function to generate the mol object for Pyscf """
        if crds is None:
            crds = np.zeros((len(atomids), 3))
        mol = gto.M(atom=[[atom, crd] for atom, crd in zip(atomids, crds)],
        basis = basis, unit='Bohr')
        return mol
