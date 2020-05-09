from collections.abc import MutableMapping
from copy import deepcopy
import os
import shutil
import subprocess
import numpy as np
#
from jinja2 import Template
#
from qctools import generate_filereader, Event
from qctools.events import join_events
#
from pysurf.colt import Colt
from pysurf.utils import exists_and_isfile
from pysurf.system import ATOMID_TO_NAME
from pysurf.logger import Logger, get_logger
#
from . import AbinitioBase
from ...system import Molecule


class ADC2questions(Colt):
    _questions = """
    ricc2_exe = ricc2 :: str
    dscf_exe = dscf :: str
    """


turbomole = """
[MP2Energy]
grep = Final MP2 energy
split = 5 :: float

[NStates]
grep = nstates=
split = -1 :: int
"""


def get_ex_energies(txt):
    ex = []
    for line in txt:
        ex += [float(line.split()[1])]
    return ex

def get_gradient(txt):
    grad = []
    for line in txt:
        line = line.replace('D', 'E')
        ls = line.split()
        grad += [[float(ls[0]), float(ls[1]), float(ls[2])]]
    return grad

def get_state(txt):
    return int(txt.split('_')[-1])

ADCExEnergies = Event('ADCExEnergies',
                    'xgrep', {'keyword': 'excitation_energies_ADC(2)',
                             'ishift' : 1,
                             'ilen' : 'nstates'},
                    func=get_ex_energies)


ExGradient = Event('ExGradient',
                    'xgrep', {'keyword': 'gradient',
                             'ishift' : 1,
                             'ilen' : 'natoms'},
                    func=get_gradient)

GsGradient = Event('GsGradient',
                    'xgrep', {'keyword': 'cycle',
                             'ishift' : 1,
                             'ilen' : 'natoms'},
                    func=get_gradient)

ExPropState = Event('ExPropState',
                    'grep', {'keyword': 'exstprop',
                             'ishift' : 0,
                             'ilen' : 1},
                    func=get_state)
#[fosc]
#grep = Strength :: 1 :: 0
#split =  2 :: float
#settings = multi=true
#
#[transmom]
#grep = Trans. Mom. :: 1 :: 0
#split = 2, 4, 6 :: float
#settings = mulit=true
#"""


#def length(kwargs):
#    natoms = kwargs['natoms']
#    if natoms % 6 == 0:
#        return (natoms // 6) * 4
#    return ((natoms // 6) + 1) * 4


#def get_gradient(scftxt):
#    gradient = []
#    for i in range(0, len(scftxt), 4):
#        x = map(float, scftxt[i + 1].split()[1:])
#        y = map(float, scftxt[i + 2].split()[1:])
#        z = map(float, scftxt[i + 3].split()[1:])
#        gradient += [[x1, y1, z1] for x1, y1, z1 in zip(x, y, z)]
#    return gradient
#




#SCFGradient = Event('SCFGradient',
#            'xgrep', {'keyword': 'Gradient of SCF Energy',
#                      'ilen': length,
#                      'ishift': 1,},
#            func=get_gradient,
#)
#
#
#CisGradient = Event('CisGradient',
#            'xgrep', {'keyword': 'Gradient of the state energy (including CIS Excitation Energy)',
#                      'ilen': length,
#                      'ishift': 1,},
#            func=get_gradient,
#)


# change by hand the events, to clear them up!
TurbomoleReader = generate_filereader("TurbomoleReader", turbomole)
TurbomoleReader.add_event("ExGradient", ExGradient)
TurbomoleReader.add_event("ExPropState", ExPropState)
grad_ex = TurbomoleReader._events['ExGradient']
prop_st = TurbomoleReader._events['ExPropState']
#tm = QChemReader._events['transmom']
grad_state = join_events(prop_st, grad_ex)
TurbomoleReader.add_event("ADCExEnergies", ADCExEnergies)
TurbomoleReader.add_event("GsGradient", GsGradient)
#QChemReader.add_event("SCFGradient", SCFGradient)
TurbomoleReader.add_event("Grad", grad_state)


tpl_energy = Template(
"""
$denconv   0.10000000E-06
$ricc2
  adc(2)
$excitations
  irrep=a mult=1 nexc={{states}} npre={{states}} nstart={{states}}
  {{fosc}}
""")

tpl_fosc = """
  spectrum  states=all operators=diplen
"""

tpl_gradient_gs = """
$denconv   0.10000000E-06
$ricc2
  geoopt model=mp2
"""

tpl_gradient_ex = Template(
"""
$denconv   0.10000000E-06
$ricc2
  adc(2)
$excitations
  irrep=a mult=1 nexc={{states}} npre={{states}} nstart={{states}}
  xgrad states=(a {{gradstates}})
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


class Turbomole(AbinitioBase):

    _questions = """
    exe_define = define :: str 
    method = ADC(2) :: str :: [ADC(2)]
    basis = cc-pVDZ
    """

    _method = {
               'ADC(2)': ADC2questions,
               }


    settings = {
        'method': 'adc2',
        'basis': 'cc-pvdz',
    }

    excited_state_settings = {
        'nstates': 3,
    }

    reader = TurbomoleReader
    #
    implemented = ['energy', 'gradient', 'fosc']

    nstates = 1

    def __init__(self, config, atomids, nstates):
        self.logger = get_logger('tm_inter.log', 'turbomole_interface')
        self.molecule = Molecule(atomids, None)
        self.nstates = nstates
        self.atomids = atomids
        self.atomnames = [ATOMID_TO_NAME[idx] for idx in atomids]
        self._update_settings(config, nstates)

            

    def _update_settings(self, config, nstates):
        self.settings = {key: config[key] for key in self.settings}
        self.excited_state_settings['nstates'] = nstates

    @classmethod
    def from_config(cls, config, atomids, nstates):
        return cls(config, atomids, nstates)

    def get(self, request):
        # update coordinates
        self.same_crd = False
        self.molecule.crd = request['crd']
        if 'energy' in request:
            self._do_energy(request)
            self.same_crd = True
        if 'gradient' in request:
            self._do_gradient(request)
        self.same_crd = False
        #
        return request

    def _do_gradient(self, request):
        grad = {}
        if 0 in request['states']:
            if exists_and_isfile('gradient'): os.remove('gradient')
            gs_gradient = self._do_gs_gradient(request)
            print(self.reader('gradient', ['GsGradient'], {'natoms': self.molecule.natoms}))
            grad[0] = np.array(self.reader('gradient', ['GsGradient'], {'natoms': self.molecule.natoms})['GsGradient'])
        fstate = True
        states = ''
        for state in request['states']:
            ex_gradient = self._do_ex_gradient(request, state)
        res = self.reader('exstates', ['ExGrad'], {'natoms': self.molecule.natoms})

        if isinstance(res['ExPropState'], list):
            for idx, state in enumerate(res['ExPropState']):
                grad[state] = np.array(res['ExGradient'][idx])
        else:
            grad[res['ExPropState']] = np.array(res['ExGradient'])
        print('Johannes grad', grad)
        request['gradient'] = grad

    def _do_energy(self, request):
        #create coord file
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        fosc = False
        mode = {'energy': ''}
        if 'fosc' in request:
            mode['fosc'] = True
        self._prepare_control(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        ricc2_output = self.submit('ricc2')
        mp2energy = self.reader(ricc2_output, ['MP2Energy'])['MP2Energy']
        nstates = self.reader('exstates', ['NStates'])['NStates']
        print('Johannes nstates', nstates)
        exenergies = self.reader('exstates', ['ADCExEnergies'], {'nstates': self.nstates-1})['ADCExEnergies']
        energies = [mp2energy]
        for ex in exenergies:
            energies += [mp2energy + ex]
        print(energies)
        request['energy'] = energies
    
    def _do_gs_gradient(self, request):
        print('Johannes do gs gradient')
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        mode = {'gradient': 0}
        self._prepare_control(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        ricc2_output = self.submit('ricc2')

    def _do_ex_gradient(self, request, state):
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        states = ''
        c=0
        mode = {'gradient': state}
        self._prepare_control(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        ricc2_output = self.submit('ricc2')

    def _prepare_control(self, mode, reffile='control_ref'):
        """write input file for qchem"""
        with open(reffile, 'r') as f:
            template = f.readlines()
        with open('control', 'w') as f:
            for line in template:
                if '$end' not in line:
                    f.write(line)
                else:
                    break
            #for energies and fosc
            if ('energy' in mode) or ('fosc' in mode):
                if 'fosc' in mode:
                    fosc = tpl_fosc
                else:
                    fosc = ''
                f.write(tpl_energy.render(states=self.nstates-1, fosc=fosc))
            #for gradients
            if 'gradient' in mode:
                #check ground state gradient calculation
                print('johannes mode grad', mode['gradient'])
                if mode['gradient'] == 0:
                    f.write(tpl_gradient_gs)
                else: 
                    f.write(tpl_gradient_ex.render(states=self.nstates-1, gradstates=mode['gradient']))
            f.write('\n$end')
    
    
    def _check_reffiles(self, coord, refcontrol=None, reforbitals=None):
        #Check whether reference control file exists
        if refcontrol is None:
            refcontrol = 'control_ref'
        if exists_and_isfile(refcontrol):
            ctrl = True
        else:
            ctrl = False
        #Check whether reference orbital file exists:
        orbs = False
        if reforbitals is None:
            if exists_and_isfile('mos'):
                orbs = True
            if exists_and_isfile('alpha') and exists_and_isfile('beta'):
                orbs = True
        else:
            if len(reforbitals) >  1:
                for reffile in reforbitals:
                    if exists_and_isfile(reffile):
                        shutil.copy(reffile, './')
                    else:
                        self.logger.error('reference orbital files not found')
            else:
                if exists_and_isfile(reffile):
                    shutil.copy(reffile, './')
                else:
                    self.logger.error('reference orbital file not found')
            orbs = True
        #run turbomole define input script if necessary
        if (orbs and ctrl) is False:
            self._run_define(coord, self.settings['basis'])

    def _write_coord(self, crd, filename='coord'):
        with open(filename, 'w') as crdfile:
            crdfile.write('$coord\n')
            for crd, atom in zip(crd, self.atomnames):
                crdfile.write(f"   {crd[0]:>19.14f}   {crd[1]:>19.14f}   {crd[2]:>19.14f}   {atom}\n")
            crdfile.write("$end")
        return filename
    
    def _run_define(self, coord, basis):
        self.logger.info('Run turbomole define script')
        define_in = f"""
        
        a {coord}
        *
        no
        b all {basis}
        *
        eht
        
        
        
        q
        """
        define = subprocess.run('define',
                                input=define_in.encode('utf-8'),
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        if not('define ended normally' in define.stderr.decode('utf-8').strip()):
            self.logger.error('Error in generating the reference control file')
        shutil.copyfile('control', 'control_ref')

    def submit(self, exe):
        output = exe + ".out"
        cmd = f"{exe} > {output}"
        os.system(cmd)
        return output
