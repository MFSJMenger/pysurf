from collections.abc import MutableMapping
from copy import deepcopy
#
from jinja2 import Template
#
from qctools import generate_filereader
#
from . import AbinitioBase


qchem = """
[SCFEnergy]
grep = Total energy in the final basis set =
split = 8 :: float

[ExcitedState]
grep = : excitation energy (eV) = :: 1 :: 1
split = 5 :: float
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


# change by hand the events, to clear them up!
QChemReader = generate_filereader("QChemReader", qchem)
QChemReader.add_event("SCFGradient", SCFGradient)
QChemReader.add_event("CisGradient", CisGradient)


tpl = Template("""
$molecule
{{chg}} {{mult}} {% for atomid, crd in mol %} 
{{mol.format(atomid, crd)}} {% endfor %}
$end

$rem {% for key, value in remsection %}
{{key}} {{value}} {% endfor %}
$end
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


class QChem(AbinitioBase):

    _questions = """
    exe = 
    remsection = :: literal
    chg = 0 :: int
    mult = 1 :: int

    """

    tpl = tpl

    reader = QChemReader

    settings = {
        'set_iter': 50,
        'exchange': 'pbe0',
        'basis': 'sto-3g', 
        'MAX_SCF_CYCLES': 500,
        'xc_grid': '000075000302',
        'mem_static': 4000,
        'mem_total': 16000,
        'SYM_IGNORE': 'true',
    }

    excited_state_settings = {
        'cis_n_roots': 3,
        'cis_singlets': True,
        'cis_triplets': False,
        'RPA': 0,
    }

    runtime = ['jobtype',]
    #
    implemented = ['energy', 'gradient']

    def __init__(self, atomids, nstates, chg, mult, qchemexe):
        self.molecule = Molecule(atomids, None)
        self.exe = qchemexe
        self.chg = chg
        self.mult = mult
        self.filename = 'qchem.in'
        self.nstates = nstates
        self.reader._events['ExcitedState'].nmax = nstates
        self._update_dct(settings)

    @classmethod
    def from_config(cls, config, atomids, nstates):
        return cls(atomids, nstates, config['chg'], config['mult'], config['exe'])

    def get(self, request):
        # update coordinates
        self.molecule.crd = request['crd']
        if 'gradient' in request:
            self._do_gradient(request)
        if 'energy' in request:
            self._do_energy(request)
        #
        return request

    def _do_gradient(self, request):
        jobtype = 'FORCE'
        remsection = deepcopy(self.settings)
        gradient = {}
        for state in request['states']:
            if state == 0:
                gradient[state] = self._do_gs_gradient(request)
            else:
                gradient[state] = self._do_ex_gradient(request, state)
        request['gradient'] = gradient

    def _do_energy(self, request):
        settings = UpdatableDict(self.settings, self.excited_state_settings)
        settings['jobtype'] = 'force'
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['SCFEnergy', 'ExcitedState'])
        request['energy'] = [out['SCFEnergy']] + out['ExcitedState']

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
        self._write_input(self.filename, settings.items())
        output = self.submit(self.filename)
        out = self.reader(output, ['CisGradient'], {'natoms': self.molecule.natoms})
        return out['CisGradient']

    def _write_input(self, filename, remsection):
        """write input file for qchem"""
        with open(filename, 'w') as f:
            f.write(self.tpl.render(chg=self.chg, mult=self.mult, 
                    mol=self.molecule, remsection=remsection))
    
    def submit(self, filename):
        output = filename + ".out"
        cmd = f"{self.exe} {filename} > {output}"
        os.system(cmd)
        return output
