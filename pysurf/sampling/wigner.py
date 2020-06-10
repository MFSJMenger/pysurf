import numpy as np
import numpy.random as random
#
from pysurf.constants import U_TO_AMU, CM_TO_HARTREE
from pysurf.molden import MoldenParser
from pysurf.spp import ModelFactory
from pysurf.colt import Colt
#
from ..system import Molecule
from ..system.atominfo import MASSES
from ..system.atominfo import ATOMNAME_TO_ID
from .normalmodes import NormalModes as nm
from .normalmodes import Mode
from .base_sampler import DynSamplerBase, DynCondition


class Molden(Colt):
    _questions = """
    moldenfile = :: existing_file
    """


class Wigner(DynSamplerBase):

    _questions = """
        # Input source for the normal modes and/or frequencies, which are used to generate the
        # initial conditions.
        from = :: str
    """

    _from = {'molden' : Molden,
             'model' : ModelFactory,
    }

    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_cases("from", {name: val.questions for name, val in cls._from.items()})
        
    
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
        elif config['from'].value == 'model':
            model = ModelFactory.instance_from_config(config['from'])
            return cls.from_model(model)
        raise Exception("only (molden, frequencies) implemented")

    @classmethod
    def from_db(cls, database):
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
    def from_molden(cls, filename):
        molden = MoldenParser(filename, ['Info', 'Freqs', 'FrCoords', 'FrNormCoords'])
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
        return cls(model, model.modes, True)

    def to_mass_weighted(self):
        if self.is_massweighted is True:
            return
        else:
            self.modes = nm.create_mass_weighted_normal_modes(self.modes, self.system)
            self.is_massweighted = True


def get_initial_condition(system, modes):
    """Generate an initial condition based on the normal modes of
       the molecule, according to a Wigner distribution
       Method is based on L. Sun, W. L. Hase J. Chem. Phys. 133, 044313 (2010).
    """
    Epot = 0.0
    #
    crd = np.copy(system.crd)
    veloc = np.zeros_like(system.crd, dtype=np.double)
    #
    for mode in modes:
        if mode.freq < 0.0:
            factor = np.sqrt(-1.0*mode.freq)
        else:
            factor = np.sqrt(mode.freq)
        #
        while True:
            # get random Q and P in the interval [-5, +5]
            Q = random.random()*10.0 - 5.0
            P = random.random()*10.0 - 5.0
            # calculate probability for this set of P and Q with Wigner distr.
            probability = wigner_gs(Q, P)
            #
            if probability > random.random():
                break  # coordinates accepted
        # now transform the dimensionless coordinate into a real one
        # QM programs directly give angular frequency (2*PI is not needed)
        # Higher frequencies give lower displacements and higher momentum.
        # Therefore scale random_Q and random_P accordingly:
        Q /= factor
        P *= factor
        # add potential energy of this mode to total potential energy
        Epot += 0.5 * mode.freq**2 * Q**2
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
    return Epot, DynSamplerBase.condition(crd, veloc, 0)


def wigner_gs(Q, P):
    """for a one-dimensional harmonic oscillator.
       Q contains the dimensionless coordinate of the
       oscillator and P contains the corresponding momentum."""
    return np.exp(-Q**2.0) * np.exp(-P**2.0)
