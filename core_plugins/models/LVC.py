import numpy as np
from pysurf.spp import Model
from pysurf.sampling.normalmodes import Mode


class LVC(Model):
    """
    Multi-state linear vibronic coupling model
    
    H = sum_j hw_j/2 * (P_j^2 + Q_j^2) + sum_ab |a> H_ab <b|
    where
    H_ab = W_ab + sum_j (kappa_j * Q_j)
    """
    
    _user_input = """
    # Config file of the LVC parameters
    parameter_file =  :: str
    """
    
    implemented = ['energy', 'gradient', 'nacs', 'Hel', 'Hel_gradient']  
    
    @classmethod
    def from_config(cls, config):
        parameter_file = config['parameter_file']
        return cls.from_parameter_file(parameter_file)
    
    @classmethod
    def from_parameter_file(cls, parameter_file):
        import configparser, json
        # Read the input
        parser = configparser.ConfigParser()
        parser.read(parameter_file)
        # Read the data
        freq = np.array(json.loads(parser['LVC']['freq']), dtype = np.float64)
        W0   = np.array(json.loads(parser['LVC']['W0']), dtype = np.float64)
        ## symmetrize
        W0 = 0.5 * (W0 + W0.T)
        #
        nstates = len(W0)
        nmodes = len(freq)
        ### 
        # Diabatic linear coupling strengths
        kappa = np.zeros((nstates, nstates, nmodes), dtype = np.float64)
        for x in parser['LVC'].keys():
            if 'kappa' in x:
                xs = x.split('_')
                i, j = int(xs[1]), int(xs[2])
                if max(i,j) > nstates - 1:
                    print(f'From the W0 matrix the number of states is {nstates}')
                kappa[i,j] = np.array(json.loads(parser['LVC'][x]))
                kappa[j,i] = kappa[i,j].copy()
        ###
        # Equilibrium value for the coordinates (optional, needed e.g. for sampling)
        Q0 = np.array(json.loads(parser['LVC']['Q0']), dtype = np.float64) \
             if 'Q0' in parser['LVC'] else np.zeros(nmodes, dtype = np.float64)
        ###
        return cls(freq, W0, kappa, Q0)
            
    def __init__(self, freq, W, kappa, Q0 = None):
        """
        Initialization of a linear vibronic coupling model

        Parameters
        ----------
        freq : array(nmodes)
            mode frequencies
        W : array(nstates,nstates)
            diabatic energies for Q = 0
        kappa : array(nstates,nstates,nmodes)
            intra- and inter-molecular coupling strengths
        Q0 : equilibrium value for the coordinates (needed, e.g., for sampling)
        """
        self.freq    = freq
        self.W       = W
        self.kappa   = kappa
        self.nstates = len(self.W)
        self.nmodes  = len(self.freq)
        
        # Attributes needed for sampling
        if Q0 is None: Q0 = np.zeros(self.nmodes, dtype = np.float64)
        self.crd    = Q0
        disp  = np.eye(self.nmodes)
        self.modes  = [Mode(w, disp[i]) for i, w in enumerate(self.freq)]
        self.masses = 1.0 / self.freq.copy()
        
    #
    # Properties
    #
    
    # Diabatic electronic Hamiltonian
    def _diab_Hel(self, Q):
        # Quadratic potential energy
        eQuad = 0.5 * np.dot(Q, (self.freq * Q))
        # Electronic Hamiltonian
        Hel = np.diag([eQuad] * self.nstates) \
            + self.W + self.kappa @ Q
        #
        return Hel
        
    # Energy
    def _energy(self, Q, adiabatic = True):
        # Electronic Hamiltonian
        Hel = self._diab_Hel(Q)
        #
        return np.linalg.eigvalsh(Hel) if adiabatic else Hel.diagonal()
    
    # Gradient of the diabatic electronic Hamiltonian
    def _diab_Hel_gradient(self, Q):
        gradient = np.zeros((self.nstates, self.nstates, self.nmodes))
        #
        wQ = self.freq * Q
        for i in range(self.nstates):
            gradient[i,i] = wQ + self.kappa[i,i]
        #
        for i in range(self.nstates):
            for j in range(i):
                gradient[i,j] = self.kappa[i,j]
                gradient[j,i] = self.kappa[j,i]
        #
        return gradient
        
    # Gradient
    def _gradient(self, Q, adiabatic = True):
        if not adiabatic:
            wQ = self.freq * Q
            return np.array([wQ + self.kappa[i,i] for i in range(self.nstates)])
        else:
            # Diabatic electronic gradient
            gradient = self._diab_Hel_gradient(Q)
            # Diabatic electronic Hamiltonian
            Hel = self._diab_Hel(Q)
            V, C = np.linalg.eigh(Hel)
            # Adiabatic gradient
            gradient_ad = np.zeros((self.nstates, self.nmodes))
            for i in range(self.nstates):
                dummy = np.zeros(self.nmodes)
                for j in range(self.nstates):
                    for k in range(self.nstates):
                        dummy += C[j,i] * C[k,i] * gradient[j,k]
                #
                gradient_ad[i] = dummy.copy()
            #
            return gradient_ad
        
    # Nonadiabatic coupling
    def _nacs(self, Q):
        # Gradient of the diabatic electronic Hamiltonian
        gradient = self._diab_Hel_gradient(Q)
        #
        Hel = self._diab_Hel(Q)
        V, C = np.linalg.eigh(Hel)
        # Nonadiabatic coupling
        nac = np.zeros((self.nstates, self.nstates, self.nmodes))
        for i in range(self.nstates):
            for j in range(self.nstates):
                dummy = np.zeros(self.nmodes)
                if i != j:
                    for k in range(self.nstates):
                        for l in range(self.nstates):
                            dummy += C[k,i] * C[l,j] * gradient[k,l]
                    dummy = - dummy / (V[i] - V[j])
                #
                nac[i,j] = dummy.copy()
        #
        return nac
                                   
    # Get the requested properties
    def get(self, request):
        """
        Parameters
        ----------
        request : 
            Request object containing the information, which properties are asked for
            
            If diabatic properties are desired (instead of the default adiabatic ones)
            the field 'adiabatic' should be given and set to False

        Returns
        -------
        request : dict
            Updated Request with the calculated properties
        """
        Q = request.crd
        adiabatic = True
        if 'adiabatic' in request.keys():
            adiabatic = request['adiabatic']
        #
        for prop in request:
            if prop == 'energy':
                energy = self._energy(Q, adiabatic = adiabatic)
                request.set('energy', energy[request.states])
            if prop == 'gradient':
                grad = self._gradient(Q, adiabatic = adiabatic)
                request.set('gradient', grad[request.states,:])
            if prop == 'nacs':
                nacs_full = self._nacs(Q)
                nacs = {}
                for i in request.states:
                    for j in request.states:
                        if i != j:
                            nacs[(i,j)] = nacs_full[i,j]
                request.set('nacs', nacs)
            if prop == 'Hel':
                request.set('Hel', self._diab_Hel(Q))
            if prop == 'Hel_gradient':
                grad_full = self._diab_Hel_gradient(Q)
                #grad = {}
                #for i in request.states:
                #    for j in request.states:
                #        grad[(i,j)] = grad_full[i,j]
                request.set('Hel_gradient', grad_full)
                #request.set('Hel_gradient', grad)
        #return None
        return request
    
    
##########################
# Test
if __name__ == "__main__":
   au2cm1 = 219474.6029    
   nmodes  = 3
   nstates = 2
   # Initialize a random LVC model
   np.random.seed(1988)
   freq  = np.random.random(nmodes) * 3500.0 / au2cm1
   W     = (np.random.random((nstates,nstates)) * 2 - 1) * 4000.0 / au2cm1
   kappa = (np.random.random((nstates,nstates,nmodes)) * 2 - 1) * 4000.0 / au2cm1
   for i in range(nstates):
       for j in range(i):
           W[i,j] = W[j,i]
           kappa[i,j] = kappa[j,i]
   # 
   lvc = LVC(freq, W, kappa)
   # or
   #lvc = LVC.from_input_file('lvc.inp')   
   #
   # TEST 1: Request energy and gradient
   Q = [0.1, -0.5, 2.1]
   from pysurf.spp.request import Request
   req = Request(Q, ['energy', 'gradient', 'nacs'], states = list(range(lvc.nstates)))   
   lvc.get(req)   
   print(req['gradient'])
   #
   # TEST 2: Sampling
   from pysurf.sampling import Wigner
   sampling = Wigner.from_model(lvc)
   print(sampling.get_condition())
   # TEST 3: Database
   from pysurf.dynamics import DynDB
   variables = DynDB._variables_model
   dimensions = {'nmodes': lvc.nmodes, 'nstates': lvc.nstates, 'nactive': lvc.nstates}
   db = DynDB.generate_database('test.db', data = variables,
                                dimensions = dimensions, model = True, sp = False)
   
    
   

        
    
