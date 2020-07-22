import numpy as np
from copy import deepcopy
from ..methodbase import Model
from pysurf.system import Mode
from pysurf.qctools.converter import energy_converter 


ev2au = 1./27.2113961

class PyrazineSala(Model):
    """ Pyrazine Models from Phys. Chem. Chem. Phys. 2014, 16, 15957 from Sala, Lasorne, Gatti and
        Guerin. The user can choose how many states are involved. Additionally the ground-state is         
        added to the models, but is not coupled. The number of modes varies with the number of states.
   """

    _questions = """
    n_states = 3 :: int :: [3, 4, 5]
    """

    q_list = ['v6a', 'v10a', 'v16a', 'v1', 'v4', 'v19b', 'v9a', 'v3', 'v8b', 'v8a', 'v5', 'v12',
            'v18a', 'v18b', 'v19a', 'v14']

    #Frequencies in cm-1 from Table 2
    w = {'v6a':  593,
         'v1':   1017,
         'v9a':  1242,
         'v8a':  1605,
         'v2':   3226,
         'v10a': 936,
         'v4':   734,
         'v5':   942,
         'v6b':  700,
         'v3':   1352,
         'v8b':  1553,
         'v7b':  3205,
         'v16a': 337,
         'v17a': 966,
         'v12':  1022,
         'v18a': 1148,
         'v19a': 1486,
         'v13':  3206,
         'v18b': 1079,
         'v14':  1364,
         'v19b': 1440,
         'v20b': 3221,
         'v16b': 417,
         'v11':  797}

    #Energies in [eV] from Table 3
    E = {0: 0.0,
         1: 3.93,
         2: 4.45,
         3: 4.79,
         4: 5.38}

    #Kappa values [eV] from Table 4
    k = {'v6a': {1: -0.081, 2: -0.168, 3: 0.128,  4:-0.184},
         'v1':  {1: -0.038, 2: -0.083, 3: -0.183, 4: -0.117},
         'v9a': {1: 0.117,  2: -0.071, 3: 0.045,  4: 0.165},
         'v8a': {1: -0.087, 2: -0.465, 3: 0.026,  4: 0.172},
         'v2':  {1: 0.022,  2: 0.060,  3: 0.018,  4: 0.030}}

    #lambda values according to Table 5 in [eV]
    l = {'v10a': 0.195,
         'v4':   0.060,
         'v5':   0.053,
         'v6b':  0.000,
         'v3':   0.065,
         'v8b':  0.219,
         'v7b':  0.020,
         'v16a': 0.112,
         'v17a': 0.018,
         'v12':  0.207,
         'v18a': 0.090,
         'v19a': 0.094,
         'v13':  0.000,
         'v18b': 0.044,
         'v14':  0.044,
         'v19b': 0.072,
         'v20b': 0.000}

    #List of allowed transitions according to the symmetry
    #Coupling state 1 and 2: B3g symmetry
    l12 = ['v6b', 'v3', 'v8b', 'v7b']
    #Coupling state 1 and 3: B1g symmetry
    l13 = ['v10a']
    #Coupling state 2 and 3: B2g symmetry
    l23 = ['v4', 'v5']
    #Coupling state 1 and 4: B1u symmetry
    l14 = ['v12', 'v18a', 'v19a', 'v13']
    #Coupling state 2 and 4: B2u symmetry
    l24 = ['v18b', 'v14', 'v19b', 'v20b']
    #Coupling state 3 and 4: Au symmetry
    l34 = ['v16a', 'v17a']

    #gamma values according to Table 5 in [eV]
    g = {'v10a': {1: -0.012, 2: -0.048, 3: -0.012, 4: -0.013},
         'v4':   {1: -0.030, 2: -0.031, 3: -0.031, 4: -0.027},
         'v5':   {1: -0.014, 2: -0.026, 3: -0.026, 4: -0.009},
         'v6b':  {1: -0.013, 2: -0.013, 3: -0.005, 4: -0.006},
         'v3':   {1: -0.006, 2: 0.006,  3: 0.001,  4: -0.004},
         'v8b':  {1: -0.012, 2: -0.012, 3: 0.007,  4: -0.043},
         'v7b':  {1: 0.003,  2: 0.003,  3: 0.004,  4: 0.003},
         'v16a': {1: 0.013,  2: -0.013, 3: -0.008, 4: -0.008},
         'v17a': {1: -0.016, 2: -0.041, 3: -0.012, 4: -0.012},
         'v12':  {1: -0.006, 2: -0.022, 3: -0.006, 4: -0.006},
         'v18a': {1: -0.006, 2: -0.002, 3: -0.005, 4: -0.006},
         'v19a': {1: -0.006, 2: -0.010, 3: -0.002, 4: -0.006},
         'v13':  {1: 0.005,  2: 0.004,  3: 0.004,  4: 0.004},
         'v18b': {1: -0.001, 2: -0.003, 3: -0.002, 4: -0.003},
         'v14':  {1: -0.019, 2: -0.021, 3: -0.020, 4: -0.021},
         'v19b': {1: -0.013, 2: -0.006, 3: -0.015, 4: -0.006},
         'v20b': {1: 0.005,  2: 0.003,  3: 0.004,  4: 0.003}}

    implemented = ["energy", "gradient", "fosc"]

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        self.nstates = config['n_states']
        if config['n_states'] == 3:
            self.q = ['v6a', 'v10a', 'v1', 'v9a', 'v8a']
            self.states = [0, 1, 3]
        if config['n_states'] == 4:
            self.q = ['v6a', 'v10a', 'v1', 'v4', 'v9a', 'v3', 'v8b', 'v8a', 'v5']
            self.states = [0, 1, 2, 3]
        if config['n_states'] == 5:
            self.q = ['v6a', 'v10a', 'v16a', 'v1', 'v4', 'v19b', 'v9a', 'v3', 'v8b',
                      'v8a', 'v5', 'v12', 'v18a', 'v18b', 'v19a', 'v14']
            self.states = [0, 1, 2, 3, 4]

        self.conv_ev = energy_converter.get_converter('eV', 'au')
        self.conv_cm = energy_converter.get_converter('cminv', 'au')

        self.frequencies = np.array([self.conv_cm(self.w[mode]) for mode in self.q])
        self.masses = np.ones(self.frequencies.size)/self.frequencies
        self.displacements = np.identity(len(self.q))
        self.modes = [Mode(freq, dis) for freq, dis in zip(self.frequencies, self.displacements)]
        self.crd = np.zeros(len(self.q))


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
        # eigh gives eigenvalues and eigenvectors in ascending order for hermitian matrix
        res = np.linalg.eigh(diab_en)
        adiab_en = res[0]
        T = res[1]
        return adiab_en, T

    def adiab_fosc(self, crd, T=None):
        """adiab_en returns a one dimensional vector with the adiabatic energies at
           the position crd
        """
        if T is None:
            diab_en = self.diab_en(crd)
            T = np.linalg.eigh(diab_en)[1]
        adiab_fosc = T.dot(self.diab_fosc(crd))
        return adiab_fosc

    def diab_en(self, crd):
        """diab_en returns the full diabatic matrix of the system
        """
        diab_en = np.zeros((self.nstates, self.nstates), dtype=float)
        #Build diabatic hamiltonian
        for state_idx, state in enumerate(self.states):
            #add excitation energy
            diab_en[state_idx, state_idx] += self.conv_ev(self.E[state])

            for idx, mode in enumerate(self.q):
                #add reference frequency terms
                diab_en[state_idx, state_idx] += self.conv_cm(self.w[mode])/2. * crd[idx]**2
       
                #Add kappa terms
                if mode in self.k and state != 0:
                    diab_en[state_idx, state_idx] += self.conv_ev(self.k[mode][state]) * crd[idx]

                #Add gamma terms
                if mode in self.g and state != 0:
                    diab_en[state_idx, state_idx] += self.conv_ev(self.g[mode][state]) * crd[idx]**2

        #Add off diagonal coupling terms:
        if self.nstates == 3:
            for idx, mode in enumerate(self.q):
                #note that the 2 state model contains the states 1 and 3!
                if mode in self.l13:
                    diab_en[1, 2] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[2, 1] = diab_en[1, 2]
        else:
            for idx, mode in enumerate(self.q):
                #note that the 2 state model contains the states 1 and 3!
                if mode in self.l13:
                    diab_en[1, 3] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[3, 1] = diab_en[1, 3]

        if self.nstates >= 4:
            for idx, mode in enumerate(self.q):
                if mode in self.l12:
                    diab_en[1, 2] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[2, 1] = diab_en[1, 2]
                if mode in self.l23:
                    diab_en[2, 3] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[3, 2] = diab_en[2, 3]
        
        if self.nstates >= 5:
            for idx, mode in enumerate(self.q):
                if mode in self.l14:
                    diab_en[1, 4] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[4, 1] = diab_en[1, 4]
                if mode in self.l24:
                    diab_en[2, 4] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[4, 2] = diab_en[2, 4]
                if mode in self.l34:
                    diab_en[3, 4] += self.conv_ev(self.l[mode]) * crd[idx]
                    diab_en[4, 3] = diab_en[3, 4]

        return diab_en

    def diab_fosc(self, crd):
        if self.nstates == 3:
            return np.array([0, 0, 1])
        if self.nstates == 4:
            return np.array([0, 0, 0, 1])
        if self.nstates == 5:
            return np.array([0, 0, 0, 1, 0])

    def adiab_grad(self, crd):
        """adiab_grad calculates and returns the numeric adiabatic gradients
        """
        adiab_grad = np.empty((self.nstates, len(self.q)), dtype=float)
        dq = 0.01
        for i in range(len(self.q)):
            crd1 = deepcopy(crd)
            crd1[i] = crd1[i] + dq
            en1, _ = self.adiab_en(crd1)
            crd2 = deepcopy(crd)
            crd2[i] = crd2[i] - dq
            en2, _ = self.adiab_en(crd2)
            adiab_grad[:, i] = (en1 - en2)/2.0/dq
        return {i: adiab_grad[i] for i in range(self.nstates)}



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
    crd = np.zeros()
    model = PyrazineSala({'':{'n_states': 3}})
    print(model.frequencies)
