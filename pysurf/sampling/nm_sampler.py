import numpy as np

from pysurf.colt import Colt
from pysurf.spp import ModelFactory
from pysurf.molden import MoldenParser
from pysurf.constants import U_TO_AMU, CM_TO_HARTREE

from .base_sampler import CrdSamplerBase
from .normalmodes import NormalModes as nm
from ..system import Molecule
from ..system.atominfo import ATOMNAME_TO_ID, MASSES
from .normalmodes import Mode
from .base_sampler import CrdCondition
from .n_grid_iter import NGridIterator

class Moldenfile(Colt):
    _questions = """
    moldenfile = :: existing_file
    """

class NMSampler(CrdSamplerBase):
    _questions = """
    # Where do the coorindates come from? If it is a moldenfile, modes are converted into dimensionless 
    # normal modes. If they come from a model, they are used as they are.
    from = molden :: str :: [molden, model]

    stepsize = 0.3 :: float

    include_combinations = False :: bool

    # Decide whether sampling should be done along all normal modes
    select_nmodes = all :: str :: [all, select]

    """

    _from = {'molden': Moldenfile,
             'model': ModelFactory,
             }

    _select_nmodes = {
            'all': "",
            'select': "mode_list = :: ilist"
            }

    _step = 0
    _mode = 0 
    _sign = 1

    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_cases("from", {name: method.questions
                                 for name, method in cls._from.items()})
        questions.generate_cases("select_nmodes", {name: value
                                 for name, value in cls._select_nmodes.items()})


    def __init__(self, config, system, modes, start=0):
        self.system = system
        self.modes = modes
        self.config = config
        if config['select_nmodes'] == 'all':
            self.sel_modes = modes
        else:
            self.sel_modes = []
            for idx in config['select_nmodes']['mode_list']:
                self.sel_modes += [modes[idx]]
        self.nmodes_sel = len(self.sel_modes)
        self.config = config
        self.stepsize = config['stepsize']
        if config['from'].value == 'model':
            self.model = True
        else:
            self.model = False
        
        self._check_modes()
        if config['include_combinations']:
            self.myiter = NGridIterator(len(self.sel_modes))

        if start != 0:
            for i in range(start):
                self.get_condition()
        

    def get_init(self):
        """Return all infos needed for the initial condition parser"""
        return {'system': self.system,
                'modes': self.modes}

    def get_condition(self):
        if self.config['include_combinations']:
            return self.get_condition_combined()
        else:
            return self.get_condition_pure()


    def get_condition_combined(self):
        vec = next(self.myiter)
        crd = np.copy(self.system.crd)

        for fac, mode in zip(vec, self.sel_modes):
            crd += np.array(fac*self.stepsize) * np.array(mode.displacements)
        print(crd)
        return CrdCondition(crd)
        

    def get_condition_pure(self):
        """Return a single created initial condition"""
        crd = np.copy(self.system.crd)
        # for reference point:
        if self._step == 0:
            self._step += 1
            return CrdCondition(crd)
        
        crd += self._sign * self._step * self.stepsize * np.array(self.sel_modes[self._mode].displacements)

        # to sample in both directions
        if self._sign == 1:
            self._sign = -1
            return CrdCondition(crd)

        if self._mode == self.nmodes_sel - 1:
            self._sign = 1
            self._mode = 0
            self._step += 1
        else:
            self._sign = 1
            self._mode += 1
        cond = CrdCondition(crd)
        return cond

    @classmethod
    def from_config(cls, config, start=0):
        """ """
        if config['from'] == 'molden':
            return cls.from_molden(config, start)
        elif config['from'].value == 'model':
            return cls.from_model(config, start)

    @classmethod
    def from_molden(cls, config, start=0):
        filename = config['from']['moldenfile']
        molden = MoldenParser(filename, ['Info', 'Freqs', 'FrCoords', 'FrNormCoords'])
        # get molecule info
        atoms = [atom for atom, _, _, _ in molden['FrCoords']]
        atomids = np.array([ATOMNAME_TO_ID[atom] for atom in atoms])
        crd = np.array([[x, y, z] for _, x, y, z in molden['FrCoords']])
        masses = np.array([MASSES[idx]*U_TO_AMU for idx in atomids])
        # create molecule
        molecule = Molecule(atomids, crd, masses)
        #
        print('nmsampler, molecule', molecule)
        modes = [Mode(freq * CM_TO_HARTREE, np.array(molden['FrNormCoords'][imode]))
                 for imode, freq in enumerate(molden['Freqs'])]
        #
        modes = nm.create_mass_weighted_normal_modes(modes, molecule)
        #
        return cls(config, molecule, modes, start=start)

    @classmethod
    def from_model(cls, config, start=0):
        model = ModelFactory.plugin_from_config(config['from']['model'])
        return cls(config, system=model, modes=model.modes, start=start)

    def _check_modes(self):
        img = [mode.freq for mode in self.modes if mode.freq < 0.0]
        nimg_freq = len(img)
        if nimg_freq == 0:
            return

        def to_strg(number):
            return "%12.8f" % number

        print(f"Found {nimg_freq} imaginary frequencies:")
        print("[" + ", ".join(map(to_strg, img)) + "]")


       

