import numpy as np

import pynof
from pysurf.spp.qm import AbinitioBase
from jinja2 import Template #build templates
from pysurf.system import Molecule

#geometry in units bohr
tpl = Template("""
units bohr
{{chg}} {{mult}} {% for atomid, crd in mol %} 
 {{mol.format(atomid, crd)}} {% endfor %}

""")

class IntPynof(AbinitioBase):
    """Interface for pynof code, which is open source
    """

    _user_input = """
    # Chage
    chg = 0 :: int
    # Multiplicity
    mult = 1 :: int
    # Basis set
    basis = cc-pvdz :: str :: cc-pvdz, 6-31G*
    # Functional
    ipnof = 8 :: int :: 4, 5, 7, 8
    # Electronic repulsion integrals
    eri = True :: bool :: True, False
    gradients = True :: bool :: True, False
    """
    tpl = tpl

    # implemented has to be overwritten by the individual classes for the methods
    implemented = ['energy', 'gradient']

    def __init__(self, config, atomids, basis, chg, mult, ipnof, eri, gradients):
        self.molecule = Molecule(atomids, None) 
        self.natoms = len(atomids)
        self.basis = basis
        self.chg = chg
        self.mult = mult
        self.ipnof = ipnof
        self.eri = eri
        self.gradients = gradients
        if not self.gradients:
            raise SystemExit("Wrong choice of gradients: it must be set to True")
        self.icall = 0
        self._last_crd = None

    @classmethod
    def from_config(cls, config, atomids):
        return cls(config, atomids, config['basis'], config['chg'], config['mult'], config['ipnof'], config['eri'], config['gradients'])


    def get(self, request):
        if self.icall == 0:
            self.read_mos = False
            self.icall = 1
        else:
            self.read_mos = True

        # Update coordinates
        self.molecule.crd = request.crd

        # Check if coordinates are the same as last call
        if self._last_crd is None or not np.allclose(self._last_crd, request.crd):
            self._do_pynof_ene_grad()
            self._last_crd = np.copy(request.crd)

        # Output requested properties
        if 'energy' in request:
            self._out_energy(request)
        if 'gradient' in request:
            self._out_gradient(request)

        return request

    def _do_pynof_ene_grad(self):
        string_geo = self.tpl.render(chg=self.chg, mult=self.mult,
                     mol=self.molecule)
        mol = pynof.molecule(string_geo)
        p = pynof.param(mol,self.basis)
        p.ipnof = self.ipnof
        p.RI = self.eri
        if self.read_mos:
            E,C,gamma,fmiug0,g = pynof.compute_energy(mol,p,C,gradients=self.gradients)
        else:
            E,C,gamma,fmiug0,g = pynof.compute_energy(mol,p,gradients=self.gradients)
        """Saving energy and gradient for the ground state"""
        self.energy = E
        self.grad = g.reshape(-1, 3)

    def _out_energy(self, request):
        """Energy of the ground state"""
        out_ene = self.energy
        request.set('energy', out_ene)

    def _out_gradient(self, request):
        """Gradientof the ground state"""
        out_grad = self.grad
        request.set('gradient', out_grad)

if __name__=='__main__':
    from pysurf.database import PySurfDB
    from pysurf.spp.request import Request
    from numpy import copy 

    db_file = "sampling.db"
    db = PySurfDB.load_database(db_file, read_only=True)
    crd = copy(db['crd'][0])
    atomids = copy(db['atomids'])

    out = IntPynof.from_questions(config="spp.inp",atomids=atomids)
    
    # Create request for ground state (state 0)
    request = Request(crd, ['energy', 'gradient'], [0])

    # Run calculation
    response = out.get(request)

    # Print results
    print("Energy:", response['energy'])
    print("Gradient:", response['gradient'][0])
