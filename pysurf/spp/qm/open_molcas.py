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
import h5py #read h5 files 
import re #regular expression syntax: string searching and manipulation

openmolcas = """
#Here is the place for grep fucntions to read the output files
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

def get_xms_energies(scftxt):
    ene = []
    for line in scftxt:
        struct = line.split()[1:]
        ene.append(float(struct[5]))
    return ene

def get_gradient(scftxt):
    grad = []
    for line in scftxt:
        struct = line.split()[1:]
        grad.append([float(struct[0]),float(struct[1]),float(struct[2])])
    return grad

def get_nacs_index(scftxt):
    nac_index = []
    for line in scftxt:
        struct = re.findall(r'\d+', line)
        nac_index = copy([int(struct[0]),int(struct[1])])
    return nac_index

def get_nacs(scftxt):
    nac = []
    for line in scftxt:
        struct = line.split()[1:]
        nac.append([float(struct[0]),float(struct[1]),float(struct[2])])
    return nac

def get_fosc(scftxt):
    ene = []
    for line in scftxt:
        struct = line.split()[1:]
        ene.append(float(struct[1]))
    return ene

XMS_Energies_OpMol = Event('XMS_Energies',
            'xgrep', {'keyword': 'XMS-CASPT2 Root',
                      'ilen': length_ene,
                      'ishift': 0,},
            func=get_xms_energies,
)

Gradient_OpMol = Event('Gradient',
            'xgrep', {'keyword': 'Molecular gradients',
                      'ilen': length,
                      'ishift': 8,},
            func=get_gradient,
)

NACs_Index_OpMol = Event('NACs_Index',
            'xgrep', {'keyword': 'Lagrangian multipliers are calculated for states no.',
                      'ilen': 1,
                      'ishift': 0,},
             settings = {"multi":True, "reset":True},
            func=get_nacs_index,
)

NACs_OpMol = Event('NACs',
            'xgrep', {'keyword': 'Total derivative coupling',
                      'ilen': length,
                      'ishift': 8,},
             settings = {"multi":True, "reset":True},
            func=get_nacs,
)

Fosc_OpMol = Event('Fosc',
            'xgrep', {'keyword': 'Dipole transition strengths',
                      'ilen': length_osc,
                      'ishift': 6,},
             settings = {"multi":True, "reset":True},
            func=get_fosc,
)
# change by hand the events, to clear them up!
OpenMolcasReader = generate_filereader("OpenMolcasReader", openmolcas)
OpenMolcasReader.add_event("Gradient_OpMol", Gradient_OpMol)
OpenMolcasReader.add_event("NACs_Index_OpMol", NACs_Index_OpMol)
OpenMolcasReader.add_event("NACs_OpMol", NACs_OpMol)
OpenMolcasReader.add_event("XMS_Energies_OpMol", XMS_Energies_OpMol)
OpenMolcasReader.add_event("Fosc_OpMol", Fosc_OpMol)


tpl = Template("""
&GATEWAY
COORD
 {{natoms}}
 Bohr
 {% for atomid, crd in mol %} 
 {{mol.format(atomid, crd)}} {% endfor %}
 {% for key, value in remsection %}
 {{key}} = {{value}} {% endfor %}

>> Do while
&SEWARD

&SCF
""")

rasscf = Template("""

&RASSCF
 {% for key, value in rasscf %}
 {{key}} = {{value}} {% endfor %}
 expert

""")


grad = Template("""

&ALASKA
 root = {{root}}
 pnew
 show
 verbose

"""
)

dc = Template("""

&ALASKA
 nac = {{nac}}
 show
 verbose

"""
)

grad_caspt2 = Template("""

&CASPT2
 {% for key, value in caspt2 %}
 {{key}} = {{value}} {% endfor %}
 rlxRoot = {{root}}

&ALASKA
"""
)

dc_caspt2 = Template("""

&CASPT2
 {% for key, value in caspt2 %}
 {{key}} = {{value}} {% endfor %}
 nac   = {{nac}}

&ALASKA
 nac = {{nac}}
 show
 verbose

"""
)

f_osc_caspt2 = Template("""

&CASPT2
 {% for key, value in caspt2 %}
 {{key}} = {{value}} {% endfor %}

 >>COPY $Project.JobMix JOB001
&RASSI
 Ejob

"""
)

f_osc_rasscf = Template("""

&RASSCF
 {% for key, value in rasscf %}
 {{key}} = {{value}} {% endfor %}
 expert

 >>COPY $Project.JobIph JOB001
&RASSI
 Ejob

"""
)

end_do = Template("""
>> EndDo
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


class OpenMolcas(AbinitioBase):

    _user_input = """
    exe = pymolcas :: str 
    basis = cc-pVDZ :: str
    group = Nosym :: str
    charge = 0 :: int
    spin = 1 :: int
    levshft = 1.0 :: float
    nactel = 4 :: int
    ras2 = 3 :: int
    frozen = 0 :: int
    deleted = 0
    ciroot = 3 3 1 :: ilist
    thrs = 1.0e-09 1.0e-05 1.0e-05 :: flist
    couplings =  yes :: str :: yes, not 
    caspt2 = :: str
    f_osc = :: str
    [caspt2(yes)]    
    xMult = all :: str
    imag = 0.20 :: float
    ipea = 0.00 :: float
    convergence = 1.0e-07 :: float
    thresholds  = 1.0e-9 1.0e-07 :: flist
    [caspt2(not)]    
    cal_caspt2 = false :: bool
    [f_osc(yes)]
    cal_f_osc = true :: bool
    [f_osc(not)]
    cal_f_osc = false :: bool
    """

    tpl = tpl

    rasscf = rasscf

    grad = grad

    dc = dc

    grad_caspt2 = grad_caspt2

    dc_caspt2 = dc_caspt2

    f_osc_caspt2 = f_osc_caspt2

    f_osc_rasscf = f_osc_rasscf

    end_do = end_do

    reader = OpenMolcasReader

    settings = {
        'group': 'Nosym',
    }

    rasscf_settings = {
        'charge': 0,
        'spin': 1,
        'levshft': 0.5,
        'nactel': 4, 
        'ras2': 3,
        'frozen': 0,
        'deleted': 0,    
        'ciroot': [3, 3, 1],
        'thrs': [1.0e-08, 1.0e-04, 1.0e-04],
    }

    xms_caspt2_settings = {
        'xMult'   : all,
        'imag'    : 0.20,
        'ipea'    : 0.00,
        'convergence' : 1.0e-07,
        'thresholds'  : [1.0e-9, 1.0e-07],
    }

    grad_settings = {
        'root': 1,
    }

    nacs_settings = {
        'nac': [1,2],
    }

    
    implemented = ['energy', 'gradient', 'nacs', 'fosc']

    def __init__(self, config, atomids, nstates, charge, spin, exe, basis, couplings, caspt2, f_osc):
        self.molecule = Molecule(atomids, None)
        self.natoms = len(atomids)
        self.exe = exe
        self.chg = charge
        self.mult = spin
        self.couplings = couplings
        self.filename = 'omolcas.in'
        self.nstates = nstates
        self._update_settings(config)
        self.caspt2 = config['caspt2']
        self.f_osc = config['f_osc']
        self.basis = basis 
        self.icall = 0
            

    def _update_settings(self, config):
        self.settings = {key: config[key] for key in self.settings}
        self.settings.update({'basis': config['basis']})
        self.rasscf_settings = {key: self.no_square_brakets(config[key]) for key in self.rasscf_settings}
        self.xms_caspt2_settings = {key: self.no_square_brakets(value) for key, value  in config['caspt2'].items()} 

    @classmethod
    def from_config(cls, config, atomids, nstates, nghost_states):
        return cls(config, atomids, nstates, config['charge'], config['spin'], config['exe'], config['basis'], config['couplings'], config['caspt2'], config['f_osc'])

    def get(self, request):
        if self.icall == 0:
            self.icall = 1
        else:
            self.rasscf_settings.update({'file': 'omolcas.RasOrb'})       
        # update coordinates
        self.molecule.crd = request.crd 
        if 'gradient' in request:
            self._out_gradient(request)
        if 'energy'  in request:
            self._out_energy(request)
        if 'nacs'  in request:
            self._out_nacs(request)
        return request

    def _read_energies(self, output, key):
        """Load keywords in memory of a hdf5 file specified by filename"""
        db =  h5py.File(output, "r")
        ene = copy(db[key])
        return ene[:self.nstates] 

    def _read_xms_energies(self, output):
        """Read energies from log file """
        out = self.reader(output, ['XMS_Energies_OpMol'], {'nstates': self.nstates})
        if not isinstance(out['XMS_Energies_OpMol'], list):
            ene = [out['XMS_Energies_OpMol']]
        else:
            ene = out['XMS_Energies_OpMol']
        return ene  

    def _read_f_osc(self, output):
        """Read osc. strengths from log file """
        out = self.reader(output, ['Fosc_OpMol'], {'nstates': self.nstates})
        if not isinstance(out['Fosc_OpMol'][0], list):
            out_f_osc = [0.] +  [out['Fosc_OpMol'][0]]
        else:
            out_f_osc = [0.] +  [value for value in out['Fosc_OpMol'][0]]
        return out_f_osc[:self.nstates] #just interactions with the groud states 

    def _read_nacs(self, output):
        out = self.reader(output, ['NACs_Index_OpMol','NACs_OpMol'], {'natoms': self.molecule.natoms})
        if not isinstance(out['NACs_OpMol'][0][0], list):
            nac = [out['NACs_OpMol']]
        else:
            nac = out['NACs_OpMol']
        nacs = {}
        leng = int(self.nstates*(self.nstates-1)/2)
        if self.nstates > 2:
            for i in range(leng):
                nacs.update({(int(out['NACs_Index_OpMol'][i][0]-1),int(out['NACs_Index_OpMol'][i][1]-1)):np.array(nac[i])})
                nacs.update({(int(out['NACs_Index_OpMol'][i][1]-1),int(out['NACs_Index_OpMol'][i][0]-1)):-np.array(nac[i])})
        elif self.nstates == 2:
            nacs.update({(int(out['NACs_Index_OpMol'][0][0]-1),int(out['NACs_Index_OpMol'][0][1]-1)):np.array(nac[0])})
            nacs.update({(int(out['NACs_Index_OpMol'][0][1]-1),int(out['NACs_Index_OpMol'][0][0]-1)):-np.array(nac[0])})
        else:
            raise SystemExit("The number of states must be higher than 1")
        return nacs

    def _read_grads(self, state):
        self.outputs = self._do_sa_casscf_ene_grad_nacs(state)
        out = self.reader(self.outputs.output_grad_nacs, ['Gradient_OpMol'], {'natoms': self.molecule.natoms})        
        if not isinstance(out['Gradient_OpMol'][0][0], list):
            grad = [out['Gradient_OpMol']]
        else:
            grad = out['Gradient_OpMol']
        return np.array(grad) 

    def _do_sa_casscf_ene_grad_nacs(self, state):
        settings = UpdatableDict(self.settings)
        rasscf_settings = UpdatableDict(self.rasscf_settings)
        root = int(state + 1)
        if self.caspt2 == "yes":
            xms_caspt2_settings = UpdatableDict(self.xms_caspt2_settings)
            self._write_caspt2_input(self.filename, settings.items(), rasscf_settings.items(), xms_caspt2_settings.items(), root)
        elif self.caspt2 == "not":
            self._write_input(self.filename, settings.items(), rasscf_settings.items(), root)
        outputs = self.submit(self.filename)
        return outputs 

    def _out_energy(self, request):
        if self.caspt2 == "yes":
            if 'fosc' in request:
                outputs = self._do_sa_casscf_ene_grad_nacs(0)
                out_ene = self._read_xms_energies(outputs.output_grad_nacs)
                out_f_osc = self._read_f_osc(outputs.output_grad_nacs) 
                request.set('energy', out_ene)
                request.set('fosc', out_f_osc)
            else:
                out_ene = self._read_xms_energies(self.outputs.output_grad_nacs)
                request.set('energy', out_ene)
        elif self.caspt2 == "not":
            if 'fosc' in request:
                outputs = self._do_sa_casscf_ene_grad_nacs(0)
                out_ene = self._read_energies(outputs.output_ene,'ROOT_ENERGIES')
                out_f_osc = self._read_f_osc(outputs.output_grad_nacs) 
                request.set('energy', out_ene)
                request.set('fosc', out_f_osc)
            else:
                out_ene = self._read_energies(self.outputs.output_ene,'ROOT_ENERGIES')
                request.set('energy', out_ene)

    def _out_gradient(self, request):
        out_gradient = {}
        for state in request.states:
            out_gradient[state] = self._read_grads(state)
        request.set('gradient', out_gradient)

    def _out_nacs(self, request):
        out_nacs = self._read_nacs(self.outputs.output_grad_nacs)
        request.set('nacs', out_nacs)
    
    def no_square_brakets(self, value):
        if isinstance(value, list):
            y = ", ".join( repr(e) for e in value)
            return y
        return value

    def _write_input(self, filename, remsection, rassection, gradsection):
        """write input file for open_molcas"""
        with open(filename, 'w') as f:
            if self.f_osc == 'yes':
                f.write(self.tpl.render(natoms=self.natoms, 
                        mol=self.molecule, remsection=remsection))
                f.write(self.f_osc_rasscf.render(rasscf=rassection))
            else:
                f.write(self.tpl.render(natoms=self.natoms, 
                        mol=self.molecule, remsection=remsection))
                f.write(self.rasscf.render(rasscf=rassection))
                f.write(self.grad.render(root=gradsection))
                if self.couplings == 'yes':
                    nac = []
                    for i in range(self.nstates):
                        for j in range(self.nstates):
                            if j>i:
                                nac = [i+1,j+1]
                                f.write(self.dc.render(nac=self.no_square_brakets(nac)))
            f.write(self.end_do.render())

    def _write_caspt2_input(self, filename, remsection, rassection, caspt2section, gradsection):
        """write input file for open_molcas"""
        with open(filename, 'w') as f:
            if self.f_osc == 'yes':
                f.write(self.tpl.render(natoms=self.natoms, 
                        mol=self.molecule, remsection=remsection))
                f.write(self.rasscf.render(rasscf=rassection))
                f.write(self.f_osc_caspt2.render(caspt2=caspt2section))
            else:
                f.write(self.tpl.render(natoms=self.natoms, 
                        mol=self.molecule, remsection=remsection))
                f.write(self.rasscf.render(rasscf=rassection))
                f.write(self.grad_caspt2.render(caspt2=caspt2section, root=gradsection))
                if self.couplings == 'yes':
                    nac = []
                    for i in range(self.nstates):
                        for j in range(self.nstates):
                            if j>i:
                                nac = [i+1,j+1]
                                f.write(self.dc_caspt2.render(caspt2=caspt2section, nac=self.no_square_brakets(nac)))
            f.write(self.end_do.render())
        
    def submit(self, filename): 
        infile = filename.replace(".in","")
        output_ene = infile + ".rasscf.h5"
        output_grad_nacs = infile + ".log"
        nt = "-f"
        cmd = f"$MOLCAS/{self.exe} {filename} {nt}"
        os.system(cmd)
        outputs = namedtuple("outputs","output_ene output_grad_nacs") 
        return outputs(output_ene, output_grad_nacs)


if __name__=='__main__':
    from pysurf.database import PySurfDB
    from pysurf.spp import SurfacePointProvider

    #def read_h5(filename, key):
    #    db =  h5py.File(filename, "r")
    #    return copy(db[key])
                      
    db_file = "init.db"
    #db_file = "sampling.db"
    db = PySurfDB.load_database(db_file, read_only=True)
    #db_file = "omolcas.rasscf.h5"
    #ene= read_h5(db_file,'CENTER_COORDINATES') 
    #print(ene)
    crd = copy(db['crd'][0])
    atomids = copy(db['atomids'])
    natoms = len(crd)
    nstates = 3
    print(crd)
    #out = OpenMolcasReader("omolcas.log",["Gradient_OpMol", "NACs_Index_OpMol", "NACs_OpMol"],{"natoms":natoms}) 
    #print("NACs_index:",out["NACs_Index_OpMol"][0][1])
    #print("NACs_index:",out["NACs_Index_OpMol"])
    #print("NACs:",out["NACs_OpMol"])
    #print("Grad:",out["Gradient_OpMol"])
    state = copy(db['currstate']) + 1
    #out = OpenMolcas.from_questions(atomids, nstates = 3, config = "spp.inp", nghost_states = None)
    #spp = SurfacePointProvider.from_questions(['energy', 'gradient', 'nacs'], nstates, natoms, atomids=atomids, config='spp.inp')
    spp = SurfacePointProvider.from_questions(['energy', 'fosc'], nstates, natoms, atomids=atomids, config='spp.inp')
    #res = spp.request(crd, ['energy','gradient','nacs'], states=[state])
    res = spp.request(crd, ['energy','fosc'], states=[1])
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
