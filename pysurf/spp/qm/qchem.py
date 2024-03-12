from collections.abc import MutableMapping
from collections import namedtuple
from copy import deepcopy
import os
import numpy as np
#
from jinja2 import Template
#
from qctools import generate_filereader, Event
from qctools.events import join_events
#
from . import AbinitioBase
from ...system import Molecule
#
from scipy.special import comb
from itertools import combinations

qchem = """
[SCFEnergy]
grep = Total energy in the final basis set =
split = 8 :: float

[ExcitedState]
grep = : excitation energy (eV) = :: 1 :: 1
split = 5 :: float
settings = multi=true

[SF_S_2]
grep = <S**2>     : :: 1 :: 0
split = 2 :: float
settings = multi=true

[fosc]
grep = Strength :: 1 :: 0
split =  2 :: float
settings = multi=true

[transmom]
grep = Trans. Mom. :: 1 :: 0
split = 2, 4, 6 :: float
settings = multi=true

[nacs]
grep = between states :: 1 :: 0
split = 2, 4 :: int
settings = multi=true, reset=True

[EOMEE]
grep = EOMEE transition :: 1 :: 1
split = 3 :: float
settings = multi=true
"""


def length(kwargs):
    natoms = kwargs['natoms']
    if natoms % 6 == 0:
        return (natoms // 6) * 4
    return ((natoms // 6) + 1) * 4


def get_gradient(scftxt):
    gradient = []
    for i in range(0, len(scftxt), 4):
        x = map(float, scftxt[i + 1].split()[1:])
        y = map(float, scftxt[i + 2].split()[1:])
        z = map(float, scftxt[i + 3].split()[1:])
        gradient += [[x1, y1, z1] for x1, y1, z1 in zip(x, y, z)]
    return gradient


def get_nacs(scftxt):
    nac = []
    for line in scftxt:
        struct = line.split()[1:]
        nac.append([float(struct[0]),float(struct[1]),float(struct[2])])
    return nac

def sf_nstates(kwargs):
    nstates = kwargs['nstates']
    return int(comb(nstates, 2)) 

def get_sf_fosc(scftxt):
    fosc = {}
    for line in scftxt:
        osc = line.split()[0:]
        fosc.update({(int(osc[0]),int(osc[1])):float(osc[5])})
    return fosc


SCFGradient = Event('SCFGradient',
            'xgrep', {'keyword': 'Gradient of SCF Energy',
                      'ilen': length,
                      'ishift': 1,},
            func=get_gradient,
)


CisGradient = Event('CisGradient',
            'xgrep', {'keyword': 'Gradient of the state energy (including CIS Excitation Energy)',
                      'ilen': length,
                      'ishift': 1,},
            func=get_gradient,
)

EOMGradient = Event('EOMGradient',
            'xgrep', {'keyword': 'Full Analytical Gradient',
                      'ilen': length,
                      'ishift': 1,},
            func=get_gradient,
)

NACouplingEOMCCSD = Event('NACouplingEOMCCSD',
            'xgrep', {'keyword': 'NAC d^x_IJ (CI part), a.u.',
                      'ilen': "natoms",
                      'ishift': 3,},
             settings = {"multi":True, "reset":True},
            func=get_nacs,
)

NACouplingGEx = Event('NACouplingGEx',
            'xgrep', {'keyword': 'DC between ground and excited states with ETF',
                      'ilen': "natoms",
                      'ishift': 3,},
             settings = {"multi":True, "reset":True},
            func=get_nacs,
)

NACouplingEx = Event('NACouplingEx',
            'xgrep', {'keyword': 'CIS derivative coupling with ETF',
                      'ilen': "natoms",
                      'ishift': 3,},
             settings = {"multi":True, "reset":True},
            func=get_nacs,
)

NACouplingSFEx = Event('NACouplingSFEx',
            'xgrep', {'keyword': 'SF-CIS derivative coupling with ETF',
                      'ilen': "natoms",
                      'ishift': 3,},
             settings = {"multi":True, "reset":True},
            func=get_nacs,
)

Sf_Fosc = Event('Sf_Fosc',
            'xgrep', {'keyword': 'Transition Moments Between Excited States',
                      'ilen': sf_nstates,
                      'ishift': 4,},
             settings = {"multi":True, "reset":True},
            func=get_sf_fosc,
)

# change by hand the events, to clear them up!
QChemReader = generate_filereader("QChemReader", qchem)
ex_st = QChemReader._events['ExcitedState']
sf_st = QChemReader._events['SF_S_2']
osc = QChemReader._events['fosc']
tm = QChemReader._events['transmom']
eom_ex_st = QChemReader._events['EOMEE']
together = join_events(ex_st, tm, osc)
sf_together = join_events(ex_st, sf_st)
QChemReader.add_event('ExcitedStateInfo', together)
QChemReader.add_event('Sf_ExcitedState', sf_together)
QChemReader.add_event('Sf_Fosc', Sf_Fosc)
QChemReader.add_event("SCFGradient", SCFGradient)
QChemReader.add_event("CisGradient", CisGradient)
QChemReader.add_event("NACouplingGEx", NACouplingGEx)
QChemReader.add_event("NACouplingEx", NACouplingEx)
QChemReader.add_event("NACouplingSFEx", NACouplingSFEx)
QChemReader.add_event("EOMGradient", EOMGradient)
QChemReader.add_event("NACouplingEOMCCSD", NACouplingEOMCCSD)


tpl = Template("""
$molecule
{{chg}} {{mult}} {% for atomid, crd in mol %} 
{{mol.format(atomid, crd)}} {% endfor %}
$end

$rem {% for key, value in remsection %}
{{key}} {{value}} {% endfor %}
$end
""")

dc = Template("""

$derivative_coupling 
0 is the reference state
{{coupled}}
$end
"""
)

bs = Template("""

{{basis_set}}
"""
)


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


class QChem(AbinitioBase):

    _user_input = """
    exe = qchem :: str 
    nthreads = 4 :: int
#    remsection = :: literal
    chg = 0 :: int
    mult = 1 :: int
    method = :: str
    basis = cc-pVDZ :: str
    input_bohr = true :: bool
    set_iter = 150 :: int
    mem_static = 4000 :: int
    mem_total = 16000 :: int
    sym_ignore = true :: bool
    couplings = nacs :: str :: nacs, wf_overlap
    spin_flip = :: str
    [method(tddft)]
    exchange = pbe0 :: str
    max_scf_cycles = 500 :: int
    xc_grid = 000075000302
    dft_d = EMPIRICAL_GRIMME :: str :: false, EMPIRICAL_GRIMME, EMPIRICAL_CHG, EMPIRICAL_GRIMME3
    scf_algorithm  = diis :: str :: diis, diis_gdm, rca_diis
    thresh_diis_switch  = 4 :: int 
    [method(eom-ccsd)]
    method = eom-ccsd :: str
    ee_triplets = 5 :: ilist 
    cc_symmetry = false :: bool   
    scf_algorithm  = diis_gdm :: str :: diis, diis_gdm, rca_diis
    cc_trans_prop = 2 :: int 
    eom_davidson_convergence = 9 :: int   
    scf_convergence = 9 :: int            
    cc_convergence = 9 :: int
    [spin_flip(true)]
    spin_flip = true :: bool
    cis_s2_thresh = 120 :: int
    sts_mom = true :: bool :: true, false
    [spin_flip(false)]
    spin_flip = false :: bool
    """

    tpl = tpl

    dc = dc

    bs = bs

    reader = QChemReader


    settings = {
        'set_iter': 150,
        'mem_static': 4000,
        'mem_total': 16000,
        'sym_ignore': 'true',
        'input_bohr': True,
    }

    excited_state_settings = {
        'cis_n_roots': 3,
        'cis_singlets': True,
        'cis_triplets': False,
        'RPA': 0,
    }

    sf_fosc_excited_state_settings = {
        'sts_mom' : True,
    }

    eomccsd_grad_settings = {
        'cc_state_to_opt': [1,3],
    }

    nonadiabatic_coupligs_settings = {
        'calc_nac': True,
    }


    wf_overlap_settings = {
       'scf_convergence' : 7,
       'scf_guess' : 'core',
       'wf_overlap_minsurf': 1,
       'wf_overlap_nsurf': 5,
       'qmmm_ext_gifs': 1,
       'dump_wf_overlap': 1,
    }


    runtime = ['jobtype',]
        
    #
    implemented = ['energy', 'gradient', 'fosc', 'transmom', 'nacs', 'wf_overlap']

    def __init__(self, config, atomids, nstates, chg, mult, exe, nthreads, couplings, spin_flip, basis, method):
        self.molecule = Molecule(atomids, None)
        self.exe = exe
        self.nthreads = nthreads
        self.chg = chg
        self.mult = mult
        self.filename = 'qchem.in'
        self.nstates = nstates
        self.reader._events['ExcitedStateInfo'].nmax = nstates -1
        self._update_settings(config)
        self.couplings = couplings 
        self.method = config['method']
        self.basis = basis 
        self.icall = 0
        if self.couplings == "wf_overlap":
            self.wf_overlap_settings['wf_overlap_nsurf'] = nstates -1
            

    def _update_settings(self, config):
        self.settings = {key: config[key] for key in self.settings}
        self.settings.update({'basis': config['basis']})
        self.settings.update({key: value for key, value in config['spin_flip'].items()})
        self.settings.update({key: value for key, value in config['method'].items()})
        self.spin_flip = config['spin_flip']['spin_flip']
        self.sts_mom = config['spin_flip']['sts_mom']
        #for key, value in config['method'].items():
        #    self.settings[key] = value

    @classmethod
    def from_config(cls, config, atomids, nstates, nghost_states):
        return cls(config, atomids, nstates, config['chg'], config['mult'], config['exe'], config['nthreads'], config['couplings'], config['spin_flip'], config['basis'], config['method'])

    def get(self, request):
        if self.icall == 0:
            self.wf_overlap_settings['dump_wf_overlap'] = 1
            self.icall = 1
        else:
            self.wf_overlap_settings['dump_wf_overlap'] = 2
        # update coordinates
        self.molecule.crd = request.crd 
        if self.method == 'tddft':
            if 'gradient' in request:
                if self.spin_flip != True:
                    self._do_gradient(request)
                else:
                    self._do_sf_gradient(request)
            if 'energy' in request:
                self._do_energy(request)
            if 'nacs' in request:
                if self.spin_flip != True:
                    self._do_nacs(request)
                else:
                    self._do_sf_nacs(request) 
            if 'wf_overlap' in request:
                self._do_energy(request)
            return request
        elif self.method == 'eom-ccsd':
            if 'gradient' in request:
                self._do_eomccsd_gradient(request)
            if 'energy' in request:
                self._do_eomccsd_energy(request)
            if 'nacs' in request:
                self._do_eomccsd_nacs(request)
            return request


    def cou_states(self):
        y = ""
        for x in range(self.nstates):
            y = y + str(x) + " "
        return y

    def sf_cou_states(self, index):
        y = ""
        for x in index:
            y = y + str(x) + " "
        return y

    def _do_gradient(self, request):
        gradient = {}
        for state in request.states:
            if state == 0:
                gradient[state] = self._do_gs_gradient(request)
            else:
                gradient[state] = self._do_ex_gradient(request, state)
        request.set('gradient', gradient)

    def _do_sf_gradient(self, request):
        gradient = {}
        for state in request.states:
            gradient[state] = self._do_sf_ex_gradient(request, state + 1)
        request.set('gradient', gradient)

    def _do_eomccsd_gradient(self, request):
        gradient = {}
        for state in request.states:
            gradient[state] = self._do_ex_eomccsd_gradient(request, state + 1)
        request.set('gradient', gradient)

    def ov_matrix(self): 
        ov_states = self.nstates - 1
        save =  open("save", "r+")
        for line in save:
            if "QCFILEPREF:" in line:
                path = line.split()[-1]
        out = path + "/398.0"
        file_id = open(out, "rb")
        w = np.fromfile(file_id);
        w = w.reshape(int(ov_states),int(ov_states))
        w = np.pad(w,(1,0)) #adding zeros for the coupling with the ground state
        #np.set_printoptions(precision=2)
        return w

    def basis_gen(self): 
        file = open("basis_gen.ini", "r")       
        for line in file:
            if "costum_basis" in line:
                w = file.read()
        file.close() 
        return w

    def _do_eomccsd_energy(self, request):
        settings = UpdatableDict(self.settings)
        settings['jobtype'] = 'sp'
        settings['ee_triplets'] = [self.nstates]
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['EOMEE'])
        if not isinstance(out['EOMEE'], list):
            outst = [out['EOMEE']]
        else:
            outst = [en[0] for en in out['EOMEE']]
        request.set('energy', outst)
        

    def _do_energy(self, request):
        if self.couplings in ("nacs","semi_coup"):
            settings = UpdatableDict(self.settings, self.excited_state_settings)
            settings['jobtype'] = 'sp'
            if self.spin_flip != True:
                settings['cis_n_roots'] = self.nstates - 1
                self._write_input(self.filename, settings.items())
                output = self.submit(self.filename)
                out = self.reader(output, ['SCFEnergy', 'ExcitedStateInfo'])
                if not isinstance(out['ExcitedState'], list):
                    outst = [out['ExcitedState']]
                else:
                    outst = [en[0] for en in out['ExcitedState']]
                request.set('energy', [out['SCFEnergy'][0]] + outst)
                if 'fosc' in request:
                    if not isinstance(out['fosc'], list):
                        outfosc = [0.] + out['fosc']
                    else:
                        outfosc = [0.] + [value[0] for value in out['fosc']]
                    request.set('fosc', outfosc)
            else:
                settings['cis_n_roots'] = self.nstates + 1
                if self.sts_mom == True:
                    settings = UpdatableDict(settings, self.sf_fosc_excited_state_settings)
                self._write_input(self.filename, settings.items())
                output = self.submit(self.filename)
                out = self.reader(output, ['Sf_ExcitedState']) 
                s2_threshold = float(settings['cis_s2_thresh']/100)
                self.index = [i for i, s in enumerate(out['SF_S_2']) if s[0] < s2_threshold] 
                sf_ene = [out['ExcitedState'][idx] for idx in self.index]
                if int(self.nstates - len(self.index)) == 0:
                    self.ntriples = 1 #one triple state
                    if not isinstance(sf_ene, list):
                        outst = [sf_ene]
                    else:
                        outst = [en[0] for en in sf_ene]
                    request.set('energy', outst)
                elif int(self.nstates - len(self.index)) < 0:
                    self.ntriples = 0 #zero triple states  
                    if not isinstance(sf_ene, list):
                        outst = [sf_ene]
                    else:
                        outst = [en[0] for en in sf_ene]
                    request.set('energy', outst[:self.nstates])
                else:
                    raise SystemExit("High spin contamination in more than one electronic states")
                    #self.nstates = self.nstates + 1
                if 'fosc' in request:
                    sf_out = self.reader(output, ['Sf_Fosc'], {'nstates': self.nstates + 1})
                    sf_combs = [i for i in list(combinations(self.index,2))[:len(sf_ene)-1]]
                    sf_combs = np.array(sf_combs)
                    sf_osc = [sf_out['Sf_Fosc'][0][i+1,j+1] for i,j in sf_combs]
                    if not isinstance(sf_osc, list):
                        outfosc = [0.] + [value for value in sf_osc]
                    else:
                        outfosc = [0.] + sf_osc
                    request.set('fosc', outfosc)
        elif self.couplings == "wf_overlap":
            settings = UpdatableDict(self.settings, self.excited_state_settings, self.wf_overlap_settings)
            settings['jobtype'] = 'sp'
            settings['cis_n_roots'] = self.nstates
            settings['wf_overlap_nsurf'] = self.nstates - 1
            self._write_input(self.filename, settings.items())
            output = self.submit_save(self.filename)
            out = self.reader(output, ['SCFEnergy', 'ExcitedStateInfo'])
            if not isinstance(out['ExcitedState'], list):
                outst = [out['ExcitedState']]
            else:
                outst = [en[0] for en in out['ExcitedState']]
            request.set('energy', [out['SCFEnergy'][0]] + outst)
        if 'transmom' in request:
            if not isinstance(out['transmom'], list):
                outtransmom = [0.] + [out['transmom']]
            else:
                outtransmom = [0.] + out['transmom']
            request.set('transmom', outtransmom)
        if 'wf_overlap' in request:    
            wf_ov = self.ov_matrix()
            request.set('wf_overlap', wf_ov)


    def _do_gs_gradient(self, request):
        settings = UpdatableDict(self.settings)
        settings['jobtype'] = 'force'
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['SCFGradient'], {'natoms': self.molecule.natoms})
        return out['SCFGradient']

    def _do_ex_gradient(self, request, state):
        settings = UpdatableDict(self.settings, self.excited_state_settings)
        settings['jobtype'] = 'force'
        settings['CIS_STATE_DERIV'] = state
        settings['cis_n_roots'] = self.nstates
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['CisGradient'], {'natoms': self.molecule.natoms})
        return out['CisGradient']

    def _do_sf_ex_gradient(self, request, state):
        settings = UpdatableDict(self.settings, self.excited_state_settings)
        settings['jobtype'] = 'force'
        if state == 1:
            settings['CIS_STATE_DERIV'] = state
            settings['cis_n_roots'] = self.nstates + 1
        elif state > 1:
            settings['CIS_STATE_DERIV'] = state + 1
            settings['cis_n_roots'] = self.nstates + 1
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['CisGradient'], {'natoms': self.molecule.natoms})
        return out['CisGradient']

    def _do_ex_eomccsd_gradient(self, request, state):
        settings = UpdatableDict(self.settings, self.eomccsd_grad_settings)
        settings['jobtype'] = 'force'
        settings['ee_triplets'] = [self.nstates]
        settings['cc_state_to_opt'] = [1,state]
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['EOMGradient'], {'natoms': self.molecule.natoms})
        return out['EOMGradient']

    def _do_nacs(self, request):
        settings = UpdatableDict(self.settings, self.excited_state_settings, self.nonadiabatic_coupligs_settings)
        settings['cis_der_numstate'] = self.nstates
        self._write_input_nacs(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['nacs','NACouplingGEx','NACouplingEx'], {'natoms': self.molecule.natoms})
        if not isinstance(out['NACouplingGEx'][0][0], list):
            nac_g = [out['NACouplingGEx']]
        else:
            nac_g = out['NACouplingGEx']
        if not isinstance(out['NACouplingEx'][0][0], list):
            nac_ex = [out['NACouplingEx']]
        else:
            nac_ex = out['NACouplingEx']
        nac_t = nac_g + nac_ex 
        nacs = {}
        leng = int(self.nstates*(self.nstates-1)/2)
        if self.nstates > 2:
            for i in range(leng):
                nacs.update({(out['nacs'][i][0],out['nacs'][i][1]):np.array(nac_t[i])})
                nacs.update({(out['nacs'][i][1],out['nacs'][i][0]):-np.array(nac_t[i])})
            request.set('nacs', nacs)
        elif self.nstates == 2:
            nacs.update({(out['nacs'][0][0],out['nacs'][0][1]):np.array(nac_t[0])})
            nacs.update({(out['nacs'][0][1],out['nacs'][0][0]):-np.array(nac_t[0])})
            request.set('nacs', nacs)
        else:
            raise SystemExit("The number of states must be higher than 1")

    def _do_sf_nacs(self, request):
        settings = UpdatableDict(self.settings, self.excited_state_settings, self.nonadiabatic_coupligs_settings)
        settings['cis_der_numstate'] = self.nstates
        settings['cis_n_roots'] = self.nstates + 1
        self._write_input_sf_nacs(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['nacs','NACouplingSFEx'], {'natoms': self.molecule.natoms})
        if not isinstance(out['NACouplingSFEx'][0][0], list):
            nac_sfex = [out['NACouplingSFEx']]
        else:
            nac_sfex = out['NACouplingSFEx']
        sf_range = []
        for i in range(self.nstates):
            for j in range(self.nstates):
                if i < j:
                    sf_range.append([i,j])
        nacs = {}
        leng = int(self.nstates*(self.nstates-1)/2)
        if self.nstates > 2:
            for i in range(leng):
                nacs.update({(sf_range[i][0],sf_range[i][1]):np.array(nac_sfex[i])})
                nacs.update({(sf_range[i][1],sf_range[i][0]):-np.array(nac_sfex[i])})
            request.set('nacs', nacs)
        elif self.nstates == 2:
            nacs.update({(sf_range[0][0],sf_range[0][1]):np.array(nac_sfex[0])})
            nacs.update({(sf_range[0][1],sf_range[0][0]):-np.array(nac_sfex[0])})
            request.set('nacs', nacs)
        else:
            raise SystemExit("The number of states must be higher than 1")

    def _do_eomccsd_nacs(self, request):
        settings = UpdatableDict(self.settings)
        settings['ee_triplets'] = [self.nstates]
        settings['calc_nac'] = 2
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['NACouplingEOMCCSD'], {'natoms': self.molecule.natoms})
        if not isinstance(out['NACouplingEOMCCSD'][0][0], list):
            nac = [out['NACouplingEOMCCSD']]
        else:
            nac = out['NACouplingEOMCCSD']
        eom_range = []
        for i in range(self.nstates):
            for j in range(self.nstates):
                if i < j:
                    eom_range.append([i,j])
        nacs = {}
        leng = int(self.nstates*(self.nstates-1)/2)
        if self.nstates > 2:
            for i in range(leng):
                nacs.update({(eom_range[i][0],eom_range[i][1]):np.array(nac[i])})
                nacs.update({(eom_range[i][1],eom_range[i][0]):-np.array(nac[i])})
            request.set('nacs', nacs)
        elif self.nstates == 2:
            nacs.update({(eom_range[0][0],eom_range[0][1]):np.array(nac[0])})
            nacs.update({(eom_range[0][1],eom_range[0][0]):-np.array(nac[0])})
            request.set('nacs', nacs)
        else:
            raise SystemExit("The number of states must be higher than 1")
    

    def _write_input(self, filename, remsection):
        """write input file for qchem"""
        with open(filename, 'w') as f:
            f.write(self.tpl.render(chg=self.chg, mult=self.mult, 
                    mol=self.molecule, remsection=remsection))
            if self.basis == 'gen':
                basis_set = self.basis_gen()
                f.write(self.bs.render(basis_set=basis_set)) 

    def _write_input_nacs(self, filename, remsection):
        """write input file for qchem when nacs are required"""
        with open(filename, 'w') as f:
            f.write(self.tpl.render(chg=self.chg, mult=self.mult, 
                    mol=self.molecule, remsection=remsection))
            coupled = self.cou_states() 
            f.write(self.dc.render(coupled=coupled)) 
            if self.basis == 'gen':
                basis_set = self.basis_gen()
                f.write(self.bs.render(basis_set=basis_set)) 

    def _write_input_sf_nacs(self, filename, remsection):
        """write input file for qchem when nacs are required"""
        with open(filename, 'w') as f:
            f.write(self.tpl.render(chg=self.chg, mult=self.mult, 
                    mol=self.molecule, remsection=remsection))
            if self.ntriples == 1:
                index = [i+1 for i in self.index]
            elif self.ntriples == 0:
                index = [i for i in range(1, (self.nstates + 1))]                
            coupled = self.sf_cou_states(index) 
            f.write(self.dc.render(coupled=coupled)) 
            if self.basis == 'gen':
                basis_set = self.basis_gen()
                f.write(self.bs.render(basis_set=basis_set)) 
    
    def submit(self, filename): 
        output = filename + ".out"
        nt = "-nt"
        cmd = f"{self.exe} {nt} {self.nthreads} {filename} > {output}"
        os.system(cmd)
        return output

    def submit_save(self, filename): 
        output = filename + ".out"
        nt = "-nt"
        save_file = filename 
        save = "-save"
        cmd = f"{self.exe} {save} {nt} {self.nthreads} {filename} {output} {save_file} > save"
        os.system(cmd)
        return output





#if __name__=='__main__':
#    atomids = ["Cl","H"]
#    molecule = Molecule(atomids, None)
#    #print(molecule.natoms)
#    #out = QChemReader("qchem.out",["nacs","NACouplingGEx","NACouplingEx"],{"natoms":2}) 
#    #out = QChemReader("sf_1_3_4_nac_HCL.out",["nacs","NACouplingSFEx"],{"natoms":2}) 
#    #out = QChemReader("qchem.out",["NACouplingGEx","NACouplingEx"],{"natoms":2}) 
#    #out = QChemReader("sf_force_HCl.out",['Sf_ExcitedState'],{'natoms': molecule.natoms}) 
#    out = QChemReader("nac_he3.out",['NACouplingEOMCCSD'],{'natoms': molecule.natoms}) 
#    #out = QChemReader("force_HCl.out",['SCFEnergy', 'ExcitedStateInfo'],{'natoms': molecule.natoms}) 
#    #out = QChemReader("force_HCl.out",["CisGradient"],{'natoms': 2}) 
#    #out = QChem.from_questions(atomids, nstates = 5, config = "test_1.inp")
#    #out.submit_save("water_8.inp")
#    #out = QChemReader("sf_force_HCl.out",['SCFEnergy', 'ExcitedStateInfo'],{'natoms': molecule.natoms}) 
#    #out = QChemReader("sf_force_HCl.out",["ExcitedState"],{'natoms': molecule.natoms}) 
#    #out = QChemReader("sf_force_HCl.out",["ExcitedStateInfo"],{'natoms': molecule.natoms}) 
#    #print(out.basis)
#    #print(out.spin_flip)
#    #print(out['ExcitedState'])
#    #print(out['SF_S_2'])
#    #print(out['nacs'])
#    #print(out['NACouplingEx'])
#    #print(out['NACouplingGEx'])
#    print(out['NACouplingEOMCCSD'])
#    #res = out.ov_matrix()
#    #print(res)
#    #print(out["nacs"],out["NACouplingGEx"][0],out["NACouplingEx"][2])
#    #print(out["NACouplingGEx"],out["NACouplingEx"])
#    #print(out["CisGradient"])
