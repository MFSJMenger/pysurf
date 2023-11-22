import numpy as np
from numpy import random
#
from colt import Colt
#
from ..constants import U_TO_AMU, CM_TO_HARTREE
from ..molden import MoldenParser, parse_molden
from ..spp import ModelFactory
#
from ..system import Molecule
from ..system.atominfo import MASSES
from ..system.atominfo import ATOMNAME_TO_ID
from .normalmodes import NormalModes as nm
from .normalmodes import Mode
from .base_sampler import DynSamplerBase
#


class Molden(Colt):
    _user_input = """
    moldenfile = :: existing_file
    """


class Wigner(DynSamplerBase):

    _user_input = """
        # Input source for the normal modes and/or frequencies, which are used to generate the
        # initial conditions.
        from = :: str
    """

    _from = {'molden' : Molden,
             'model' : ModelFactory,
    }

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases("from", {name: val.colt_user_input
                                          for name, val in cls._from.items()})


    def __init__(self, system, modes, check=True, is_massweighted=False):
        """Initialize a new Wigner Sampling with a molecule class
           and a normal Mode class"""
        self.system = system
        self.modes = modes
        self.is_massweighted = is_massweighted
        if check is True:
            self._check_modes()

    @classmethod
    def from_config(cls, config, start=None):
        """ """
        #start keyword is not needed here, but has to be provided for DynSamplerBase
        if config['from'] == 'molden':
            return cls.from_molden(config['from']['moldenfile'])
        if config['from'].value == 'model':
            model = ModelFactory.plugin_from_config(config['from']['model'])
            return cls.from_model(model)
        raise Exception("only (molden, frequencies) implemented")

    @classmethod
    def from_db(cls, database):
        # CHECK: Is this intentional?
        return cls(None, None, check=False)

    def get_init(self):
        """Return all infos needed for the initial condition parser"""
        return {'system': self.system,
                'modes': self.modes}

    def get_condition(self):
        """Return a single created initial condition"""
        _, conds = get_initial_condition(self.system, self.modes)
        return conds

    def _check_modes(self):
        img = [mode.freq for mode in self.modes if mode.freq < 0.0]
        nimg_freq = len(img)
        if nimg_freq == 0:
            return

        def to_strg(number):
            return "%12.8f" % number

        print(f"Found {nimg_freq} imaginary frequencies:")
        print("[" + ", ".join(map(to_strg, img)) + "]")

    @classmethod
    def from_molden(cls, filename, format=2):
        if format == 1:
            molden = MoldenParser(filename, ['Info', 'Freqs', 'FrCoords', 'FrNormCoords'])
        else:
            freqs, frcoords, frnorm = parse_molden(filename)
            molden = {'FrCoords': frcoords,
                      'Freqs': freqs,
                      'FrNormCoords': frnorm,}
        # get molecule info
        atoms = [atom for atom, _, _, _ in molden['FrCoords']]
        atomids = np.array([ATOMNAME_TO_ID[atom] for atom in atoms])
        crd = np.array([[x, y, z] for _, x, y, z in molden['FrCoords']])
        masses = np.array([MASSES[idx]*U_TO_AMU for idx in atomids])
        # create molecule
        molecule = Molecule(atomids, crd, masses)
        #
        modes = [Mode(freq * CM_TO_HARTREE, np.array(molden['FrNormCoords'][imode]))
                 for imode, freq in enumerate(molden['Freqs'])]
        #
        modes = nm.create_mass_weighted_normal_modes(modes, molecule)
        #
        return cls(molecule, modes, True)

    # method to create initial conditions for model systems like the pyrazine model
    @classmethod
    def from_model(cls, model):
        """Create Wigner class from a analytic model"""
        return cls(model, model.modes, True)

    def to_mass_weighted(self):
        """transform normalmodes to massweighted"""
        if self.is_massweighted is True:
            return
        self.modes = nm.create_mass_weighted_normal_modes(self.modes, self.system)
        self.is_massweighted = True


def get_random(*args):
    if len(args) == 0:
        return random.random()
    if len(args) == 2:
        return _get_random(*args)
    raise ValueError("Please give upper and lower bounds")


def _get_random(lower, upper):
    assert upper > lower
    diff = upper - lower
    return random.random() * diff + lower


def get_initial_condition(system, modes):
    """Wigner sampling condition according to

       L. Sun, W. L. Hase J. Chem. Phys. 133, 044313 (2010).

       parameters taken from SHARC in accordance to github.com/sharc-md
       especially the bounds [-5, +5]
    """
    epot = 0.0
    #
    crd = np.copy(system.crd)
    veloc = np.zeros_like(system.crd, dtype=np.double)
    #
    for mode in modes:
        # Factor is sqrt(angular freq)
        if mode.freq < 0.0:
            factor = np.sqrt(-1.0*mode.freq)
        else:
            factor = np.sqrt(mode.freq)
        #
        while True:
            # get random Q and P in the interval [-5, +5]
            Q = get_random(-5, 5)
            P = get_random(-5, 5)
            #
            probability = wigner_gs(Q, P)
            #
            if probability > get_random():
                break  # coordinates accepted

        Q /= factor
        P *= factor
        #
        epot += 0.5 * mode.freq**2 * Q**2
        # scaling for crd, veloc sampling
        scale = np.copy(mode.displacements)
        for i, mass in enumerate(system.masses):
            scale[i] *= 1.0/np.sqrt(mass)
        #
        crd += Q * scale
        veloc += P * scale
    #
    # remove translational/rotational dofs
    #
    return epot, DynSamplerBase.condition(crd, veloc, 0)


def wigner_gs(Q, P):
    """for a one-dimensional harmonic oscillator.
       Q contains the dimensionless coordinate of the
       oscillator and P contains the corresponding momentum."""
    return np.exp(-Q**2.0) * np.exp(-P**2.0)
