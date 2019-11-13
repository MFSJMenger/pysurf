import os
#
from copy import deepcopy
#
from collections import namedtuple
#
import numpy as np
import numpy.random as random
#
from pysurf.molden import MoldenParser
#
from .initialconditions import InitialConditions



class WignerSampling(object):
    questions = """ 
        # Input source for the normal modes and/or frequencies, which are used to generate the 
        # initial conditions.
        # Possible options are:
        # - molden
        # - frequencies
        from = none :: str :: [molden, frequencies]
        
        # If initial conditions are generated from a molden file, this subquestion asks for the 
        # molden file.
        [from(molden)]
        moldenfile = none

        # If initial conditions are generated from frequencies, here a list of frequencies has to
        # be given for the modes. The list can be in the python list format or just comma separated
        # values.
        [from(frequencies)]
        frequencies = none
    """
    def __init__(self, molecule, modes, is_massweighted=False):
        """Initialize a new Wigner Sampling with a molecule class
           and a normal Mode class"""
        self.molecule = molecule
        self.modes = modes
        self.is_massweighted = is_massweighted
        self._check_modes()

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
    def from_config(cls, config):
        """ """    
        if config['from'] == 'molden':
            sampling = cls.from_molden(config['from']['filename'])
        elif config['from'] == 'freqs':
            sampling = cls.from_freqs(config['from']['freqs'])
        else:
            raise Exception("only (molden, freqs) implemented")

    @classmethod
    def from_molden(cls, filename):
        molden = MoldenParser(filename, ['Info', 'Freqs', 'FrCoords', 'FrNormCoords'])
        # get molecule info
        atoms = [atom for atom, _, _, _ in molden['FrCoords']]
        atomids = np.array([atomname_to_id[atom] for atom in atoms])
        crd = np.array([[x, y, z] for _, x, y, z in molden['FrCoords']])
        masses = np.array([MASSES[idx]*U_TO_AMU for idx in atomids])
        # create molecule
        molecule = Molecule(atomids, crd, masses)
        #
        modes = [Mode(freq * CM_TO_HARTREE, np.array(molden['FrNormCoords'][imode]))
                 for imode, freq in enumerate(molden['Freqs'])]
        #
        modes = create_mass_weighted_normal_modes(modes, molecule)
        #
        return cls(molecule, modes, True)

    # method to create initial conditions for model systems like the pyrazine model
    @classmethod
    def from_freqs(cls, freqs):
        nfreqs = len(freqs)
        masses = 1./freqs
        mol = Molecule(np.ones(nfreqs), np.zeros(nfreqs, dtype=np.double), np.ones(nfreqs, dtype=np.double))
        mol.masses = masses
        displacement = np.ones(nfreqs, dtype=np.double)
        displacement /= np.linalg.norm(displacement)
        modes = [Mode(freq, displacement) for freq in freqs]
        return cls(mol, modes, True)
        

    def create_initial_conditions(self, filename, nconditions, E_equil=0.0, model=False):
        return InitialConditions.from_conditions(filename, self.molecule, self.modes, nconditions, E_equil, model)

    def to_mass_weighted(self):
        if self.is_massweighted is True:
            return
        else:
            self.modes = create_mass_weighted_normal_modes(self.modes, self.molecule)
            self.is_massweighted = True



def get_initial_condition(molecule, modes):
    """Generate an initial condition based on the normal modes of
       the molecule, according to a Wigner distribution
       Method is based on L. Sun, W. L. Hase J. Chem. Phys. 133, 044313 (2010).
    """
    Epot = 0.0
    #
    crd = np.copy(molecule.crd)
    veloc = np.zeros((molecule.natoms, 3), dtype=np.double)
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
        for i, mass in enumerate(molecule.masses):
            scale[i] *= 1.0/np.sqrt(mass)
        #
        crd += Q * scale
        veloc += P * scale
    #
    # remove translational/rotational dofs
    #
    return Epot, InitialCondition(crd, veloc)

def compute_ekin(veloc, masses):
    e = 0.0
    for i, vel in enumerate(veloc):
        e += np.sum(0.5*(vel**2)*masses[i])
    return e


def wigner_gs(Q, P):
    """for a one-dimensional harmonic oscillator.
       Q contains the dimensionless coordinate of the
       oscillator and P contains the corresponding momentum."""
    return np.exp(-Q**2.0) * np.exp(-P**2.0)


if __name__ == '__main__':
    random.seed(16661)
    sampling = WignerSampling.from_molden('molden.dat')
    sampling.create_initial_conditions('init.db', 3, E_equil=0.0)
