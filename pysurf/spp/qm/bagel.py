from collections.abc import MutableMapping
from collections import namedtuple
import os
import numpy as np
#
from jinja2 import Template #build templates
#
from qctools import generate_filereader, Event
from qctools.events import join_events
#
from . import AbinitioBase
from ...system import Molecule
#
from scipy.special import comb
from itertools import combinations
from numpy import copy 
import re #regular expression syntax: string searching and manipulation

bagel = """
#Here is the place for grep fucntions to read the output files
[fosc]
grep = * Oscillator strength : :: 1 :: 0
split =  4 :: float
settings = multi=true
"""

def length(kwargs):
    natoms = kwargs['natoms']
    return int(natoms) 

def length_ene(kwargs):
    nstates = kwargs['nstates']
    return int(nstates) 

def length_osc(kwargs):
    nstates = kwargs['nstates']
    return int(comb(nstates,2)) 

def get_energies(scftxt):
    ene = []
    for line in scftxt:
        struct = line.split()
        ene.append(float(struct[0]))
    return ene

def get_gradient(scftxt):
    grad = []
    for line in scftxt:
        struct = line.split()[1:]
        grad.append([float(struct[0]),float(struct[1]),float(struct[2])])
    return grad

def get_nacs(scftxt):
    nac = []
    for line in scftxt:
        struct = line.split()[1:]
        nac.append([float(struct[0]),float(struct[1]),float(struct[2])])
    return nac

Energies_BAGEL = Event('Energies',
            'xgrep', {'keyword': '',
                      'ilen': length_ene,
                      'ishift': 0,},
            func=get_energies,
)

Gradient_BAGEL = Event('Gradient',
            'xgrep', {'keyword': '',
                      'ilen': length,
                      'ishift': 1,},
            func=get_gradient,
)

NACs_BAGEL = Event('NACs',
            'xgrep', {'keyword': '',
                      'ilen': length,
                      'ishift': 1,},
            func=get_nacs,
)

BAGELReader = generate_filereader("BAGELReader", bagel)
osc = BAGELReader._events["fosc"]
BAGELReader.add_event("Fosc_BAGEL", osc)
BAGELReader.add_event("Gradient_BAGEL", Gradient_BAGEL)
BAGELReader.add_event("NACs_BAGEL", NACs_BAGEL)
BAGELReader.add_event("Energies_BAGEL", Energies_BAGEL)


bagel_ini = Template("""{ "bagel" : [

""")

bagel_geo = Template("""
{
  "title" : "molecule",{% for key, value in geometry %}
  "{{key}}" : "{{value}}",{% endfor %}
  "geometry" : [{% for atomid, crd in mol %}{% if not loop.last %}
                {{mol.bagel_xyz(atomid, crd)}},{% else %}
                {{mol.bagel_xyz(atomid, crd)}} {% endif %}{% endfor %}
  ]
 },

""")


bagel_mo_load = Template("""
  {
    "title" : "load_ref",
    "file" : "guess_orbitals",
    "continue_geom" : false
  },
""")

bagel_mo_save = Template("""
  {
    "title" : "save_ref",
    "file" : "guess_orbitals"
  }
""")

bagel_osc_ini =  Template("""
  {
   "title" : "force",
   "target": 0,
   "dipole": true,
""")

bagel_grad_ini =  Template("""
  {
   "title" : "forces",
   "grads" : [
      {{state}}
""")

bagel_nacs = Template("""
      {{nacs}}
""")

bagel_end_bracket_w_coma = Template(""" 
],
""")

bagel_end_w_coma = Template(""" 
},
""")

bagel_end_wo_coma = Template(""" 
}
""")

bagel_grad_end_casscf = Template("""
   "export" : true,
   "method" : [
  {
   "title" : "casscf",{% for key, value in casscf %}{% if not loop.last %}
   "{{key}}" : {{value}},{% else %}
   "{{key}}" : {{value}} {% endif %}{% endfor %}
  }
]
""")

bagel_mo_casscf =  Template("""
  {
   "title" : "casscf",{% for key, value in casscf %}{% if not loop.last %}
   "{{key}}" : {{value}},{% else %}
   "{{key}}" : {{value}} {% endif %}{% endfor %}
  },
""")

bagel_grad_end_caspt2 =  Template("""
   "export" : true,
   "method" : [
         {
          "title" : "caspt2",
          "smith" : {
                "method" : "caspt2",{% for key, value in caspt2 %}{% if not loop.last %}
                "{{key}}" : {{value}},{% else %}
                "{{key}}" : {{value}} {% endif %}{% endfor %}
            },{% for key, value in casscf %}{% if not loop.last %}
          "{{key}}" : {{value}},{% else %}
          "{{key}}" : {{value}} {% endif %}{% endfor %}
        }
       ]
""")

bagel_end = Template("""
]}
""")

class UpdatableDict(MutableMapping):

    def __init__(self, *dicts):
        self._dicts = dicts
        self._data = {}
        self._own = []

    def clean(self):
        self._data = {}
        self._own = []

    def __setitem__(self, key, value):
        if key not in self:
            self._own.append(key)
        self._data[key] = value

    def __getitem__(self, key):
        value = self._data.get(key, None)
        if value is not None:
            return value
        for dct in self._dicts:
            value = dct.get(key, None)
            if value is not None:
                return value
        raise KeyError

    def __contains__(self, key):
        if key not in self._data:
            return any(key in dct for dct in self._dicts)
        return True

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        for key in self._own:
            yield key
        for dct in self._dicts:
            for key in dct:
                yield key

    def __len__(self):
        length = 0
        for dct in self._dicts:
            length += len(dct)
        length += len(self._own)
        return length


class Bagel(AbinitioBase):

    _user_input = """
    exe = BAGEL :: str 
    basis = cc-pvdz :: str
    df_basis = cc-pvdz-jkfit :: str
    charge = 0 :: int
    nspin = 1 :: int
    nopen = 0 :: int
    nact = 3 :: int
    nstate = 2 :: int
    nclosed = 6 :: int
    couplings =  yes :: str :: yes, not 
    active_space = :: str
    caspt2 = :: str
    f_osc = :: str  
    [caspt2(yes)]    
    ms = true :: str :: true, false
    xms = true :: str :: true, false
    sssr = true :: str :: true, false
    shift = 0.2 :: float
    frozen =  true :: str :: true, false 
    [active_space(yes)]
    active = 1 2 :: ilist
    [active_space(not)]
    active = false :: bool
    [caspt2(not)]    
    cal_caspt2 = false :: bool
    [f_osc(yes)]
    cal_f_osc = true :: bool
    [f_osc(not)]
    cal_f_osc = false :: bool
    """

    bagel_ini = bagel_ini 
    bagel_geo = bagel_geo 
    bagel_mo_load = bagel_mo_load 
    bagel_mo_save = bagel_mo_save 
    bagel_mo_casscf = bagel_mo_casscf
    bagel_grad_ini = bagel_grad_ini
    bagel_nacs = bagel_nacs
    bagel_grad_end_casscf = bagel_grad_end_casscf
    bagel_grad_end_caspt2 = bagel_grad_end_caspt2
    bagel_osc_ini = bagel_osc_ini
    bagel_end_bracket_w_coma = bagel_end_bracket_w_coma
    bagel_end_w_coma = bagel_end_w_coma
    bagel_end_wo_coma = bagel_end_wo_coma
    bagel_end = bagel_end


    reader = BAGELReader

    settings = {
        'basis' : 'cc-pvdz',
        'df_basis' : 'cc-pvdz-jkfit',
    }

    casscf_settings = {
        "nopen" : 0,
        "nact" : 3,
        "nstate" : 2,
        "nclosed" : 6,
    }

    active_space_settings = {
        "active" : [1,2]
    }

    xms_caspt2_settings = {
        "method" : "caspt2",
        "ms" : "true",
        "xms" : "true",
        "sssr" : "true",
        "shift" : 0.2,
        "frozen" : "true",
    }

    
    implemented = ['energy', 'gradient', 'nacs', 'fosc']

    def __init__(self, config, atomids, nstates, charge, nspin, exe, basis, df_basis, couplings, caspt2, f_osc, active_space):
        self.molecule = Molecule(atomids, None)
        self.natoms = len(atomids)
        self.exe = exe
        self.chg = charge
        self.mult = nspin
        self.couplings = couplings
        self.filename = 'bagel.json'
        self.nstates = nstates
        self.reader._events['Fosc_BAGEL'].nmax = nstates
        self._update_settings(config)
        self.caspt2 = config['caspt2']
        self.f_osc = config['f_osc']
        self.active = config['active_space']
        self.basis = basis 
        self.df_basis = df_basis 
        self.icall = 0
            

    def _update_settings(self, config):
        self.settings = {key: config[key] for key in self.settings}
        self.settings.update({'basis': config['basis']})
        self.settings.update({'df_basis': config['df_basis']})
        self.settings.update({'angstrom': False })
        self.casscf_settings = {key: config[key] for key in self.casscf_settings}
        self.xms_caspt2_settings = {key: value for key, value  in config['caspt2'].items()} 
        self.active_space_settings = {key: value for key, value  in config['active_space'].items()} 

    @classmethod
    def from_config(cls, config, atomids, nstates, nghost_states):
        return cls(config, atomids, nstates, config['charge'], config['nspin'], config['exe'], config['basis'], config['df_basis'], config['couplings'], config['caspt2'], config['f_osc'], config['active_space'])

    def get(self, request):
        if self.icall == 0:
            self.read_guess_mo = "not"
            self.icall = 1
        else:
            self.read_guess_mo = "yes"
        # update coordinates
        self.molecule.crd = request.crd 
        if 'gradient' in request:
            self._out_gradient(request)
        if 'energy'  in request:
            self._out_energy(request)
        if 'nacs'  in request:
            self._out_nacs(request)
        return request

    def _read_energies(self):
        """Read energies from output file """
        output = "ENERGY.out"
        out = self.reader(output, ['Energies_BAGEL'], {'nstates': self.nstates})
        if not isinstance(out['Energies_BAGEL'], list):
            ene = [out['Energies_BAGEL']]
        else:
            ene = out['Energies_BAGEL']
        return ene  

    def _read_f_osc(self, output):
        """Read osc. strengths from log file """
        out = self.reader(output, ['Fosc_BAGEL'])
        out = [value[0] for value in out['Fosc_BAGEL']]
        if not isinstance(out, list):
            out_f_osc = [0.] +  [out]
        else:
            out_f_osc = [0.] +  out
        return out_f_osc[:self.nstates] #just interactions with the groud states 

    def _nac_result(self,output):
        out = self.reader(output, ['NACs_BAGEL'], {'natoms': self.molecule.natoms})
        if not isinstance(out['NACs_BAGEL'][0][0], list):
            nac = [out['NACs_BAGEL']]
        else:
            nac = out['NACs_BAGEL']
        return nac
        
    def _read_nacs(self):
        nacs = {}
        prefixed = [filename for filename in os.listdir('.') if filename.startswith("NACME")]
        leng = int(self.nstates*(self.nstates-1)/2)
        if len(prefixed)!=leng: 
            raise SystemExit("The number of states is wrong. It must be higher than 1")
        for nac_r in prefixed:
            index = [int(s) for s in re.findall(r'\d+', nac_r)]
            nac = self._nac_result(nac_r)
            nacs.update({(index[0],index[1]):np.array(nac[0])})
            nacs.update({(index[1],index[0]):-np.array(nac[0])})
        return nacs

    def _read_grads(self, state):
        self._do_sa_casscf_ene_grad_nacs(state)
        prefixed = [filename for filename in os.listdir('.') if filename.startswith("FORCE")]
        output = prefixed[0]
        out = self.reader(output, ['Gradient_BAGEL'], {'natoms': self.molecule.natoms})        
        if not isinstance(out['Gradient_BAGEL'][0][0], list):
            grad = [out['Gradient_BAGEL']]
        else:
            grad = out['Gradient_BAGEL']
        return np.array(grad) 

    def _do_sa_casscf_ene_grad_nacs(self, state):
        settings = UpdatableDict(self.settings)
        if self.active == "yes":
            casscf_settings = UpdatableDict(self.casscf_settings, self.active_space_settings)
        else:
            casscf_settings = UpdatableDict(self.casscf_settings)
        root = int(state)
        if self.caspt2 == "yes":
            xms_caspt2_settings = UpdatableDict(self.xms_caspt2_settings)
            self._write_caspt2_input(self.filename, settings.items(), casscf_settings.items(), xms_caspt2_settings.items(), root)
        elif self.caspt2 == "not":
            self._write_input(self.filename, settings.items(), casscf_settings.items(), root)
        self.submit(self.filename)

    def _out_energy(self, request):
        if 'fosc' in request:
            self._do_sa_casscf_ene_grad_nacs(0)
            output = "bagel.out" 
            out_ene = self._read_energies()
            out_f_osc = self._read_f_osc(output) 
            request.set('energy', out_ene)
            request.set('fosc', out_f_osc)
        else:
            out_ene = self._read_energies()
            request.set('energy', out_ene)

    def _out_gradient(self, request):
        out_gradient = {}
        for state in request.states:
            out_gradient[state] = self._read_grads(state)
        request.set('gradient', out_gradient)

    def _out_nacs(self, request):
        out_nacs = self._read_nacs()
        request.set('nacs', out_nacs)
    
    def _grad_yes_no_coupling(self,root):
        if self.couplings == 'yes':
            return "{ \"title\" : \"force\", \"target\" : %s }," % root
        else:
            return "{ \"title\" : \"force\", \"target\" : %s }" % root

    def _nacs_no_last(self, nac_i, nac_j):
        return "{ \"title\" : \"nacme\", \"target\" : %(i)s, \"target2\" : %(j)s }," % {"i":nac_i, "j":nac_j}

    def _nacs_last(self, nac_i, nac_j):
        return "{ \"title\" : \"nacme\", \"target\" : %(i)s, \"target2\" : %(j)s }" % {"i":nac_i, "j":nac_j}

    def _write_input(self, filename, remsection, cassection, gradsection):
        """write input file for BAGEL"""
        with open(filename, 'w') as f:
            if self.f_osc == 'yes':
                f.write(self.bagel_ini.render()) 
                f.write(self.bagel_geo.render(geometry=remsection,mol=self.molecule))
                f.write(self.bagel_osc_ini.render()) 
                f.write(self.bagel_grad_end_casscf.render(casscf=cassection)) 
                f.write(self.bagel_end_wo_coma.render()) 
            else:
                f.write(self.bagel_ini.render()) 
                f.write(self.bagel_geo.render(geometry=remsection,mol=self.molecule))
                if self.read_guess_mo == "yes":
                    f.write(self.bagel_mo_load.render())
                f.write(self.bagel_grad_ini.render(state=self._grad_yes_no_coupling(gradsection))) 
                if self.couplings == 'yes':
                    for i in range(self.nstates):
                        for j in range(self.nstates):
                            if j>i:
                                if i == int(self.nstates-2) and j == int(self.nstates-1):
                                    f.write(self.bagel_nacs.render(nacs=self._nacs_last(i,j)))
                                else:
                                    f.write(self.dc.render(nacs=self._nacs_no_last(i,j)))
                f.write(self.bagel_end_bracket_w_coma.render()) 
                f.write(self.bagel_grad_end_casscf.render(casscf=cassection)) 
                f.write(self.bagel_end_w_coma.render()) 
                f.write(self.bagel_mo_casscf.render(casscf=cassection))
                f.write(self.bagel_mo_save.render()) 
            f.write(self.bagel_end.render())

    def _write_caspt2_input(self, filename, remsection, cassection, caspt2section, gradsection):
        """write input file for BAGEL"""
        with open(filename, 'w') as f:
            if self.f_osc == 'yes':
                f.write(self.bagel_ini.render()) 
                f.write(self.bagel_geo.render(geometry=remsection,mol=self.molecule))
                f.write(self.bagel_osc_ini.render()) 
                f.write(self.bagel_grad_end_caspt2.render(caspt2=caspt2section,casscf=cassection)) 
                f.write(self.bagel_end_wo_coma.render()) 
            else:
                f.write(self.bagel_ini.render()) 
                f.write(self.bagel_geo.render(geometry=remsection,mol=self.molecule))
                if self.read_guess_mo == "yes":
                    f.write(self.bagel_mo_load.render())
                f.write(self.bagel_grad_ini.render(state=self._grad_yes_no_coupling(gradsection))) 
                if self.couplings == 'yes':
                    for i in range(self.nstates):
                        for j in range(self.nstates):
                            if j>i:
                                if i == int(self.nstates-2) and j == int(self.nstates-1):
                                    f.write(self.bagel_nacs.render(nacs=self._nacs_last(i,j)))
                                else:
                                    f.write(self.bagel_nacs.render(nacs=self._nacs_no_last(i,j)))
                f.write(self.bagel_end_bracket_w_coma.render()) 
                f.write(self.bagel_grad_end_caspt2.render(caspt2=caspt2section,casscf=cassection)) 
                f.write(self.bagel_end_w_coma.render()) 
                f.write(self.bagel_mo_casscf.render(casscf=cassection))
                f.write(self.bagel_mo_save.render()) 
            f.write(self.bagel_end.render())
        
    def submit(self, filename): 
        infile = filename.replace(".json","")
        output = infile + ".out"
        mpi = "mpirun"
        cmd = f"{mpi} ${self.exe} {filename} > {output}"
        os.system(cmd)


if __name__=='__main__':
    from pysurf.database import PySurfDB
    from pysurf.spp import SurfacePointProvider

    #db_file = "init.db"
    db_file = "sampling.db"
    db = PySurfDB.load_database(db_file, read_only=True)
    #db_file = "omolcas.rasscf.h5"
    #ene= read_h5(db_file,'CENTER_COORDINATES') 
    #print(ene)
    crd = copy(db['crd'][0])
    atomids = copy(db['atomids'])
    natoms = len(crd)
    nstates = 3 
    state = 1
    #print(crd)
    #out = BAGELReader("omolcas.log",["Gradient_BAGEL", "NACs_Index_OpMol", "NACs_BAGEL"],{"natoms":natoms}) 
    #print("NACs_index:",out["NACs_Index_OpMol"][0][1])
    #print("NACs_index:",out["NACs_Index_OpMol"])
    #print("NACs:",out["NACs_BAGEL"])
    #print("Grad:",out["Gradient_BAGEL"])
    #state = copy(db['currstate'])
    #out = Bagel.from_questions(atomids, nstates = 2, config = "test.inp", nghost_states = None)
    #out._do_sa_casscf_ene_grad_nacs(state)
    #out._read_grads(state)
    #out._read_energies()
    #spp = SurfacePointProvider.from_questions(['energy', 'gradient', 'nacs'], nstates, natoms, atomids=atomids, config='test.inp')
    spp = SurfacePointProvider.from_questions(['energy', 'fosc'], nstates, natoms, atomids=atomids, config='test.inp')
    #res = spp.request(crd, ['energy','gradient','nacs'], states=[state])
    res = spp.request(crd, ['energy','fosc'], states=[state])
    print("ENE:",res['energy'])
    print("OSC:",res['fosc'])
    #print("NACS:",res['nacs'])
    #print("GRAD:",res['gradient'][1])

    #print(molecule.natoms)
   #out = QChemReader("qchem.out",["nacs","NACouplingGEx","NACouplingEx"],{"natoms":2}) 
    #out = QChemReader("sf_1_3_4_nac_HCL.out",["nacs","NACouplingSFEx"],{"natoms":2}) 
    #out = QChemReader("qchem.out",["NACouplingGEx","NACouplingEx"],{"natoms":2}) 
    #out = QChemReader("sf_force_HCl.out",['Sf_ExcitedState'],{'natoms': molecule.natoms}) 
    #out = QChemReader("nac_he3.out",['NACouplingEOMCCSD'],{'natoms': molecule.natoms}) 
    #out = QChemReader("force_HCl.out",['SCFEnergy', 'ExcitedStateInfo'],{'natoms': molecule.natoms}) 
    #out = QChemReader("force_HCl.out",["CisGradient"],{'natoms': 2}) 
   #result = self.spp.request(crd, ['gradient'], states=[curr_state])
   #out._write_input("hola.inp")
   #out.get("sampling.db")
    #out.submit_save("water_8.inp")
    #out = QChemReader("sf_force_HCl.out",['SCFEnergy', 'ExcitedStateInfo'],{'natoms': molecule.natoms}) 
    #out = QChemReader("sf_force_HCl.out",["ExcitedState"],{'natoms': molecule.natoms}) 
    #out = QChemReader("sf_force_HCl.out",["ExcitedStateInfo"],{'natoms': molecule.natoms}) 
    #print(out.basis)
    #print(out.spin_flip)
    #print(out['ExcitedState'])
    #print(out['SF_S_2'])
    #print(out['nacs'])
    #print(out['NACouplingEx'])
    #print(out['NACouplingGEx'])
    #print(out['NACouplingEOMCCSD'])
    #res = out.ov_matrix()
    #print(res)
    #print(out["nacs"],out["NACouplingGEx"][0],out["NACouplingEx"][2])
    #print(out["NACouplingGEx"],out["NACouplingEx"])
    #print(out["CisGradient"])
