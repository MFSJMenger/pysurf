import numpy as np

from saoovqe import SAOOVQE
from pysurf.spp.qm import AbinitioBase
from jinja2 import Template #build templates
from pysurf.system import Molecule

#geometry in units bohr
tpl = Template("""
units    bohr
 {% for atomid, crd in mol %} 
 {{mol.format(atomid, crd)}} {% endfor %}
symmetry c1
nocom
noreorient 
""")

class INTSAOOVQE(AbinitioBase):
    """ Interface for the SAOOVQE code, which is free available 
        The communication with the SPP is the get function with the request object.
        The actual calculations are performed in separate classes. The classes need to have
        a function for each property with the name do_prop, e.g. do_energy for energy calculations.
        All properties which are implemented have to be stored in the implemented property of the 
        corresponding classes
    """

    _user_input = """
    # Basis set
    basis = cc-pvdz :: str :: cc-pvdz, 6-31G*
    # Number of active electrons in the Active-Space
    nelec_active = 4 :: int 
    frozen_indices = 6 :: int 
    active_indices = 6, 9 :: ilist
    virtual_indices = 9, 43 :: ilist 
    do_oo_process = True :: str
    """
    
    tpl = tpl

    # implemented has to be overwritten by the individual classes for the methods
    implemented = ['energy', 'gradient', 'nacs']


    def __init__(self, config, atomids, nstates, basis, nelec_active, frozen_indices, active_indices, virtual_indices, do_oo_process):
        self.molecule = Molecule(atomids, None) 
        self.natoms = len(atomids)
        self.nstates = nstates
        self.basis = basis
        self.nelec_active = nelec_active
        self.frozen_indices = [i for i in range(frozen_indices)]
        self.active_indices = [i for i in range(active_indices[0], active_indices[1])] 
        self.virtual_indices = [i for i in range(virtual_indices[0], virtual_indices[1])]
        self.num_qubits = 2 * len(active_indices)
        self.do_oo_process = do_oo_process
        self.n_mo_optimized = self.virtual_indices[-1] + 1
        self.w_a = 0.5
        self.w_b = 0.5
        self.delta = 1e-5
        self.tell_me = True
        self.ucc_ansatz = ["fermionic_SAAD", "fast"][1]
        self.bohr = 0.5291772105638411 
        self.icall = 0

        
    @classmethod
    def from_config(cls, config, atomids, nstates, nghost_states):
        return cls(config, atomids, nstates, config['basis'], config['nelec_active'], config['frozen_indices'], config['active_indices'], config['virtual_indices'], config['do_oo_process'])

    
    def get(self, request):
        if self.icall == 0:
            self.read_mos = False
            self.icall = 1
        else:
            self.read_mos = True
        # update coordinates
        self.molecule.crd = request.crd
        if 'gradient' in request:
            self._out_gradient(request)
        if 'energy'  in request:
            self._out_energy(request)
        if 'nacs'  in request:
            self._out_nacs(request)
        return request

    def _do_saoovqe_ene_grad_nacs(self, state):
        string_geo = self.tpl.render(mol=self.molecule)
        saoovqe_class = SAOOVQE(string_geo,
                          self.basis,
                          self.nelec_active,
                          self.frozen_indices,
                          self.active_indices,
                          self.virtual_indices,
                          tell_me=self.tell_me,
                          w_a=self.w_a,
                          w_b=self.w_b,
                          delta=self.delta,
                          print_timings=False, # Use this if you want to compute all the timings...
                          do_oo_process=self.do_oo_process,
                          ucc_ansatz=self.ucc_ansatz)
        """Read molecular orbitals after first timestep"""
        if self.read_mos:
            saoovqe_class.c_rhf = np.loadtxt("mos_save")
        else:
            pass
        e_0, e_1 = saoovqe_class.vqe_kernel()
        """Save molecular orbitals"""
        mos_save = saoovqe_class.c_optimized
        np.savetxt("mos_save", mos_save)
        grad = []
        nac = []
        for i in range(self.natoms):
            dx,dy,dz = saoovqe_class.get_gradient(i,state)
            nx,ny,nz = saoovqe_class.get_nac(i)
            grad.append([dx*self.bohr,dy*self.bohr,dz*self.bohr])
            nac.append([-nx*self.bohr,-ny*self.bohr,-nz*self.bohr])
        self.energies = [e_0, e_1]
        self.grads = grad
        self.nacs = nac

    def _read_grads(self,state):
        self._do_saoovqe_ene_grad_nacs(state)        
        return np.array(self.grads)

    def _read_nacs(self):
        nacs = {}
        leng = int(self.nstates*(self.nstates-1)/2)
        if self.nstates == 2:
            nacs.update({(0,1):np.array(self.nacs)})
            nacs.update({(1,0):-np.array(self.nacs)})
        else:
            raise SystemExit("The number of states is different than 2")
        return nacs

    def _out_energy(self, request):
        out_ene = self.energies
        request.set('energy', out_ene)

    def _out_gradient(self, request):
        out_gradient = {}
        for state in request.states:
            out_gradient[state] = self._read_grads(state)
        request.set('gradient', out_gradient)

    def _out_nacs(self, request):
        out_nacs = self._read_nacs()
        request.set('nacs', out_nacs)


if __name__=='__main__':
    from pysurf.database import PySurfDB
    from pysurf.spp import SurfacePointProvider
    from numpy import copy 

    db_file = "sampling.db"
    db = PySurfDB.load_database(db_file, read_only=True)
    crd = copy(db['crd'][0])
    atomids = copy(db['atomids'])
    natoms = len(crd)
    nstates = 2
    state = 1

    spp = SurfacePointProvider.from_questions(['energy', 'gradient', 'nacs'], nstates, natoms, atomids=atomids, config='spp.inp')
    res = spp.request(crd, ['energy','gradient','nacs'], states=[state])
    print("ENE:",res['energy'])
    print("NACS:",res['nacs'])
    print("GRAD:",res['gradient'][1])
