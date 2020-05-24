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
    """

class DFTquestions(Colt):
    _questions = """
    functional = :: str
    grid = m4 :: str
    """


turbomole = """
[MP2Energy]
grep = Final MP2 energy
split = 5 :: float

[DSCFEnergy]
grep = |  total energy
split = 4 :: float

[ESCFEnergy]
grep = Total energy
split = 2 :: float
settings = multi=True

[ESCFFosc]
grep = Oscillator strength :: 1 :: 2
split = 2 :: float
settings = multi=True

[ADCFosc]
grep = oscillator strength
split = 5 :: float
settings = multi=True

[NStates]
grep = nstates=
split = -1 :: int
"""

# Additional functions for the reader
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

def get_fosc(txt):
    txt = txt.replace('D', 'E')
    txt_sp = txt.split()
    res = 0
    for i in [1, 3, 5]:
        res += float(txt_sp[i].strip('='))
    return res

# Events for the reader
ADCExEnergies = Event('ADCExEnergies',
                    'xgrep', {'keyword': 'excitation_energies_ADC(2)',
                             'ishift' : 1,
                             'ilen' : 'nstates'},
                    func=get_ex_energies)


ADCExGradient = Event('ADCExGradient',
                    'xgrep', {'keyword': 'gradient',
                             'ishift' : 1,
                             'ilen' : 'natoms'},
                    func=get_gradient)

Gradient = Event('Gradient',
                    'xgrep', {'keyword': 'cycle',
                             'ishift' : 'shift',
                             'ilen' : 'natoms'},
                    func=get_gradient)

ADCExPropState = Event('ADCExPropState',
                    'grep', {'keyword': 'exstprop',
                             'ishift' : 0,
                             'ilen' : 1},
                    func=get_state)


# change by hand the events, to clear them up!
TurbomoleReader = generate_filereader("TurbomoleReader", turbomole)
TurbomoleReader.add_event("ADCExGradient", ADCExGradient)
TurbomoleReader.add_event("ADCExPropState", ADCExPropState)
grad_ex = TurbomoleReader._events['ADCExGradient']
prop_st = TurbomoleReader._events['ADCExPropState']
grad_state = join_events(prop_st, grad_ex)
TurbomoleReader.add_event("ADCExEnergies", ADCExEnergies)
TurbomoleReader.add_event("Gradient", Gradient)
TurbomoleReader.add_event("ADCExGrad", grad_state)


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

tpl_dft = Template("""
$dft
  functional {{functional}}
  gridsize {{grid}}
$soes
  a {{states}}
$denconv 1d-7
$scfinstab rpas
""")


class Turbomole(AbinitioBase):

    _questions = """
    # Methode die angewendet werden soll
    method = ADC(2) :: str :: [ADC(2), DFT/TDDFT]

    basis = cc-pVDZ
    max_scf_cycles = 50 :: int
    """

    _method = {
               'ADC(2)': ADC2questions,
               'DFT/TDDFT': DFTquestions
               }


    settings = {
        'method': 'ADC(2)',
        'basis': 'cc-pvdz',
        'max_scf_cycles': 50,
    }


    reader = TurbomoleReader
    #
    implemented = ['energy', 'gradient', 'fosc']

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions for name, method in cls._method.items()})
        
    def __init__(self, config, atomids, nstates):
        self.logger = get_logger('tm_inter.log', 'turbomole_interface')
        self.molecule = Molecule(atomids, None)
        self.nstates = nstates
        self.atomids = atomids
        self.atomnames = [ATOMID_TO_NAME[idx] for idx in atomids]
        self._update_settings(config, nstates)
        self.reader._events['ESCFEnergy'].nmax = nstates

    def _update_settings(self, config, nstates):
        self.settings = {key: config[key] for key in self.settings}
        self.settings.update({key: config['method'][key] for key in config['method']})
        self.settings['nstates'] = nstates

    @classmethod
    def from_config(cls, config, atomids, nstates):
        return cls(config, atomids, nstates)

    def get(self, request):
        # update coordinates
        self.same_crd = False
        self.molecule.crd = request.crd
        if 'energy' in request:
            if self.settings['method'] == 'ADC(2)':
                self._do_energy_adc(request)
                self.same_crd = True
            elif self.settings['method'] == 'DFT/TDDFT':
                self._do_energy_dft(request)
                self.same_crd = True
        if 'gradient' in request:
            self._do_gradient(request)
        self.same_crd = False
        #
        return request

    def _do_gradient(self, request):
        grad = {}
        if 0 in request.states:
            if exists_and_isfile('gradient'): os.remove('gradient')
            if self.settings['method'] == 'ADC(2)':
                grad.update(self._do_gs_gradient_adc(request))
            elif self.settings['method'] == 'DFT/TDDFT':
                grad.update(self._do_gs_gradient_dft(request))
        for state in request.states:
            if state == 0:
                continue
            if self.settings['method'] == 'ADC(2)':
                grad.update(self._do_ex_gradient_adc(request, state))
            if self.settings['method'] == 'DFT/TDDFT':
                grad.update(self._do_ex_gradient_dft(request, state))
        request.set('gradient', grad)

    def _do_energy_adc(self, request):
        #create coord file
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        mode = {'energy': ''}
        if 'fosc' in request:
            mode['fosc'] = True
        self._prepare_control_adc(mode=mode)
        #start calculations
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        ricc2_output = self.submit('ricc2')
        #read energies
        mp2energy = self.reader(ricc2_output, ['MP2Energy'])['MP2Energy']
        energies = [mp2energy]
        if self.nstates > 1:
            exenergies = self.reader('exstates', ['ADCExEnergies'], {'nstates': self.nstates-1})['ADCExEnergies']
            if isinstance(exenergies, list):
                for ex in exenergies:
                    energies += [mp2energy + ex]
            if isinstance(exenergies, float):
                energies += [mp2energy + exenergies]
        request.set('energy', energies)
        #read oscillator strengths if requested
        if 'fosc' in request:
            fosc = [0.]
            foscres = self.reader(ricc2_output, ['ADCFosc'])['ADCFosc']
            if isinstance(foscres, list):
                fosc += foscres
            elif isinstance(foscres, float):
                fosc += [foscres]
            request.set('fosc', fosc)
    
    def _do_energy_dft(self, request):
        #create coord file
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        mode = {'energy': ''}
        self._prepare_control_dft(mode=mode)
        #start calculations
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        if self.nstates > 1:
            escf_output = self.submit('escf')
        #read energies
        dscfenergy = self.reader(dscf_output, ['DSCFEnergy'])['DSCFEnergy']
        energies = [dscfenergy]
        if self.nstates > 1:
            exenergies = self.reader(escf_output, ['ESCFEnergy'])['ESCFEnergy']
            energies = exenergies
        request.set('energy', energies)
        #read oscillator strengths if requested
        if 'fosc' in request:
            fosc = [0.]
            if self.nstates > 1:
                foscread = self.reader(escf_output, ['ESCFFosc'])['ESCFFosc']
                if isinstance(foscread, list):
                        fosc += foscread
                if isinstance(foscread, float):
                    fosc += [foscread]
            request.set('fosc', fosc)

    def _do_ex_gradient_dft(self, request, state):
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        mode = {'gradient': state}
        self._prepare_control_dft(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        grad_output = self.submit('egrad')
        return {state: np.array(self.reader('gradient', ['Gradient'], {'natoms': self.molecule.natoms, 'shift': self.molecule.natoms + 1})['Gradient'])}
    
    def _do_gs_gradient_adc(self, request):
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        mode = {'gradient': 0}
        self._prepare_control_adc(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        ricc2_output = self.submit('ricc2')
        return {0: np.array(self.reader('gradient', ['Gradient'], {'natoms': self.molecule.natoms, 'shift': self.molecule.natoms + 1})['Gradient'])}

    def _do_gs_gradient_dft(self, request):
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        mode = {'gradient': 0}
        self._prepare_control_dft(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        grad_output = self.submit('grad')
        return {0: np.array(self.reader('gradient', ['Gradient'], {'natoms': self.molecule.natoms, 'shift': self.molecule.natoms + 1})['Gradient'])}

    def _do_ex_gradient_adc(self, request, state):
        coord = self._write_coord(self.molecule.crd, filename='coord')
        self._check_reffiles(coord)
        states = ''
        c=0
        mode = {'gradient': state}
        self._prepare_control_adc(mode=mode)
        if self.same_crd is False:
            dscf_output = self.submit('dscf')
        ricc2_output = self.submit('ricc2')
        #read gradient output
        res = self.reader('exstates', ['ADCExGrad'], {'natoms': self.molecule.natoms})
        grad={}
        if isinstance(res['ADCExPropState'], list):
            for idx, state in enumerate(res['ADCExPropState']):
                grad[state] = np.array(res['ADCExGradient'][idx])
        else:
            grad[res['ADCExPropState']] = np.array(res['ADCExGradient'])
        return grad

    def _prepare_control_adc(self, mode, reffile='control_ref'):
        with open(reffile, 'r') as f:
            template = f.readlines()
        with open('control', 'w') as f:
            for line in template:
                if '$end' not in line:
                    if '$scfiterlimit' in line:
                        f.write(f"$scfiterlimit {self.settings['max_scf_cycles']}\n")
                        continue
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
                if mode['gradient'] == 0:
                    f.write(tpl_gradient_gs)
                else: 
                    f.write(tpl_gradient_ex.render(states=self.nstates-1, gradstates=mode['gradient']))
            f.write('\n$end')
    
    def _prepare_control_dft(self, mode, reffile='control_ref'):
        with open(reffile, 'r') as f:
            template = f.readlines()
        with open('control', 'w') as f:
            for line in template:
                if '$end' not in line:
                    if '$scfiterlimit' in line:
                        f.write(f"$scfiterlimit {self.settings['max_scf_cycles']}\n")
                        continue
                    #Use Turbomole default for scfdamp in DFT calculations
                    if '$scfdamp' in line:
                        f.write("$scfdamp start=0.7 step=0.050 min=0.050\n")
                        continue
                    f.write(line)
                else:
                    break
            f.write(tpl_dft.render(states=self.nstates-1, functional=self.settings['functional'], grid=self.settings['grid']))
            #for gradients
            if 'gradient' in mode:
                #check ground state gradient calculation
                if mode['gradient'] == 0:
                    pass
                else: 
                    f.write(f"$exopt {mode['gradient']}")
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
        cmd = f"{exe}"
        with open(output, 'w') as outfile:
            run = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE)
        if not('ended normally' in run.stderr.decode('utf-8')):
            self.logger.error(f"Error in Turbomole calculation: {cmd} with error: {run.stderr.decode('utf-8')}")
        return output
