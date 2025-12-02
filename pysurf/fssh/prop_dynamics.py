import numpy as np

from ..database import PySurfDB
from colt import Colt

class State(Colt):

    _user_input = """ 
    # chosen parameters
    db_file = :: existing_file 
    t = 0.0 :: float
    dt = 1.0 :: float
    mdsteps = 40000 :: float
    substeps = :: str 
    # Nose-Hoover thermostat
    thermostat = :: str
    # instate is the initial state: 0 = G.S, 1 = E_1, ...
    instate = 1 :: int
    nstates = 2 :: int
    states = 0 1 :: ilist
    ncoeff = 0.0 1.0 :: flist
    # diagonal probability is not working yet
    prob = tully :: str :: tully, lz, lz_nacs    
    rescale_vel = :: str 
    rev_vel_no_hop = true :: bool :: true, false
    coupling = nacs :: str :: nacs, wf_overlap, non_coup, semi_coup 
    method = Surface_Hopping :: str :: Surface_Hopping, Born_Oppenheimer  
    decoherence = EDC :: str :: EDC, IDC_A, IDC_S, No_DC 
    [substeps(true)]
    n_substeps = 10 :: int
    [substeps(false)]
    n_substeps = false :: bool
    [rescale_vel(momentum)]
    number_vdf = false :: str :: false, nonlinear, linear
    [rescale_vel(nacs)]
    res_nacs = true :: bool
    [thermostat(true)]
    # friction coefficient
    xi = 0.0 :: float
    # target tempertaure in Kelvin
    T = 300 :: float    
    # degrees of freedom 
    dof = nonlinear :: str :: nonlinear, linear
    [thermostat(false)]
    therm = false :: bool
    """

    def __init__(self, config, crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, rev_vel_no_hop, coupling, method, decoherence, atomids, substeps):

        self.crd = crd
        self.natoms = len(crd)
        self.atomids = atomids
        self.vel = vel
        self.mass = mass
        if model == 1:
            self.model = True
        else:
            self.model = False
        self.t = t
        self.dt = dt
        self.mdsteps = mdsteps
        self.instate = instate
        self.nstates = nstates
        self.states = states
        self.ncoeff = ncoeff
        self.prob = prob
        self.rescale_vel = config['rescale_vel'].value
        if config['rescale_vel'] == "momentum":
            self.reduced_kene = config['rescale_vel']['number_vdf']
        self.coupling = coupling
        if config['rescale_vel'] == "nacs":
            if self.coupling in ("wf_overlap, non_coup"):
                raise SystemExit("Wrong coupling method or wrong rescaling velocity approach")
        self.rev_vel_no_hop = rev_vel_no_hop
        self.method = method
        self.decoherence = decoherence
        if config['substeps'] == "true":
            self.substeps = True 
            self.n_substeps = config['substeps']['n_substeps']
        else:
            self.substeps = False 
        self.e_curr = None
        self.e_prev_step = None
        self.e_two_prev_steps = None
        self.ekin = 0
        self.epot = 0
        self.nac = {}
        self.ene = []
        self.vk = []
        self.u = []
        self.rho = []
        if np.isscalar(self.mass):
            self.natoms = 1
        elif isinstance(self.mass, np.ndarray) != True:
            self.natoms = np.array([self.mass])
        if config['thermostat'] == "true":
            self.thermostat = True 
            self.xi = config['thermostat']['xi']
            self.dof = config['thermostat']['dof']
            self.t_target = config['thermostat']['T']*3.166811e-6
            if self.dof == "nonlinear":
                self.q_eff = (3*self.natoms-6)*self.t_target*(10*self.dt)**2
            else:
                self.q_eff = (3*self.natoms-5)*self.t_target*(10*self.dt)**2
        else:
            self.thermostat = False 

    @classmethod
    def from_config(cls, config):
        crd, vel, mass, atomids, model = cls.read_db(config["db_file"])
        t = config['t']
        dt = config['dt']
        mdsteps = config['mdsteps']
        instate = config['instate']
        nstates = config['nstates']
        states = config['states']
        ncoeff = config['ncoeff']
        prob = config['prob']
        rescale_vel = config['rescale_vel']
        rev_vel_no_hop = config['rev_vel_no_hop']
        coupling = config['coupling']
        method = config['method']
        decoherence = config['decoherence']
        substeps = config['substeps']
        return cls(config, crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, rev_vel_no_hop, coupling, method, decoherence, atomids, substeps)  


    @staticmethod
    def read_db(db_file):
        db = PySurfDB.load_database(db_file, read_only=True)
        crd = np.copy(db['crd'][0])
        vel = np.copy(db['veloc'][0])
        atomids = np.copy(db['atomids'])
        mass = np.copy(db['masses'])
        model = np.copy(db['model'])
        if model == 1:
            model = True
        else:
            model = False
        return crd, vel, mass, atomids, model

    @classmethod
    def from_initial(cls, config, crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, rev_vel_no_hop, coupling, method, decoherence, atomids, substeps):
        return cls(config, crd, vel, mass, model, t, dt, mdsteps, instate, nstates, states, ncoeff, prob, rescale_vel, rev_vel_no_hop, coupling, method, decoherence, atomids, substeps)


if __name__=="__main__":
    State.from_questions(config = "prop.inp")
