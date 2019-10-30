import os
#
from copy import deepcopy
#
from collections import namedtuple
#
import numpy as np
import numpy.random as random
#
from .molden import MoldenParser
#
from .database.database import Database
from .database.dbtools import DBVariable
from .database.dbtools import DatabaseTools
from .atominfo import masses as MASSES
from .atominfo import atomname_to_id
from .constants import U_TO_AMU, CM_TO_HARTREE
#
Mode = namedtuple("Mode", ["freq", "displacements"])
#  used to save just InitialConditions (coordinates, velocities) for a given molecule
InitialCondition = namedtuple("InitialCondition", ['crd', 'veloc'])


class Molecule(object):
    """Store info of a molecule"""
    def __init__(self, atomids, crd, masses, name=None):
        self.atomids = atomids
        self.crd = crd
        self.masses = masses
        self.name = name
    @property
    def natoms(self):
        return len(self.atomids)


class WignerSampling(object):

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
        print('Johannes in from freqs')
        masses = 1./freqs
        mol = Molecule(np.ones(nfreqs), np.zeros(nfreqs, dtype=np.double), np.ones(nfreqs, dtype=np.double))
        mol.masses = masses
        displacement = np.ones(nfreqs, dtype=np.double)
        displacement /= np.linalg.norm(displacement)
        modes = [Mode(freq, displacement) for freq in freqs]
        print('Johannes modes:', modes) 
        return cls(mol, modes, True)
        

    def create_initial_conditions(self, filename, nconditions, E_equil=0.0, model=False):
        print('Johannes: model ', model)
        return InitialConditions.from_conditions(filename, self.molecule, self.modes, nconditions, E_equil, model)

    def to_mass_weighted(self):
        if self.is_massweighted is True:
            return
        else:
            self.modes = create_mass_weighted_normal_modes(self.modes, self.molecule)
            self.is_massweighted = True


class InitialConditions(object):


    def __init__(self, db, nconditions):        
        self.nconditions = nconditions
        self._db = db

    @classmethod
    def from_conditions(cls, filename, molecule, modes, nconditions=1000, E_equil=0.0, model=False):
        db = cls.create_initial_conditions(filename, molecule, modes, nconditions, E_equil, model)
        return cls(db, nconditions)

    @classmethod
    def save_condition(cls, filename, init_cond):
        natoms = len(init_cond.crd)
        nmodes = natoms*3-6
        db = Database(filename, cls.generate_settings(natoms=natoms, nmodes=nmodes, model=False))
        db.append('veloc', init_cond.veloc)
        db.append('crd', init_cond.crd)

    @classmethod
    def from_db(cls, filename, E_quil=0.0):
        db_settings = DatabaseTools.get_variables(filename, ["natoms", "nmodes"])
        db = Database(filename, cls.generate_settings(db_settings['natoms'], 
                                                      db_settings['nmodes']))
        return cls(db, db['crd'].shape[0])

    @staticmethod
    def generate_settings(natoms=0, nmodes=0, model=False):
        print('Johannes nmodes:' , nmodes)
        if model is False:
            return {
                    'dimensions': {
                        'frame': 'unlimited',
                        'nmodes': nmodes,
                        'natoms': natoms,
                        'three': 3,
                        'one': 1,
                        },
                    'variables': {
                        'modes':  DBVariable(np.double, ('nmodes', 'natoms', 'three')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'atomids': DBVariable(np.integer, ('natoms',)),
                        'masses': DBVariable(np.double, ('natoms',)),
                        'crd': DBVariable(np.double, ('frame', 'natoms', 'three')),
                        'veloc': DBVariable(np.double, ('frame', 'natoms', 'three')),
                        'energy': DBVariable(np.double, ('frame', 'three')),  # E_tot, E_pot, E_kin
                        'state': DBVariable(np.double, ('frame', 'one'))
                    },
            }

        if model is True:
            return {
                    'dimensions': {
                        'frame': 'unlimited',
                        'nmodes': nmodes,
                        'three': 3,
                        'one': 1,
                        },
                    'variables': {
                        'modes':  DBVariable(np.double, ('nmodes','nmodes')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'masses': DBVariable(np.double, ('nmodes',)),
                        'crd': DBVariable(np.double, ('frame', 'nmodes')),
                        'veloc': DBVariable(np.double, ('frame', 'nmodes')),
                        'energy': DBVariable(np.double, ('frame', 'three')),  # E_tot, E_pot, E_kin
                        'state': DBVariable(int, ('frame', 'one'))
                    },
            }

    def __iter__(self):
        self._start = 1 # skip equilibrium structure
        return self

    def get_condition(self, idx): 
        if idx >= self.nconditions:
            return None

        crd = self._db.get('crd', idx)
        veloc = self._db.get('veloc', idx)
        return InitialCondition(crd, veloc)

    @property
    def equilibrium(self):
        return self.get_condition(0)

    def __next__(self):
        cond = self.get_condition(self._start)
        if cond is not None:
            self._start += 1
            return cond
        raise StopIteration

    @property
    def molecule(self):
        return Molecule(self._db['atomids'], self._db.get('crd', 0), self._db['masses'])

    @classmethod
    def create_initial_conditions(cls, filename, molecule, modes, nconditions, E_equil, model=False):
        """ """
        if model is False:
            db = create_initial_conditions(cls.generate_settings(molecule.natoms, len(modes)),
                                       filename, molecule, modes, E_equil, model)
        else:
            db = create_initial_conditions(cls.generate_settings(nmodes=len(modes), model=True),
                                       filename, molecule, modes, E_equil, model)

        for _ in range(nconditions):
            e_pot, conds = get_initial_condition(molecule, modes)
            db.append('crd', conds.crd)
            db.append('veloc', conds.veloc)
            e_kin = compute_ekin(conds.veloc, molecule.masses)
            db.append('energy', np.array([e_kin+e_pot, e_pot, e_kin]))
            db.increase
        return db


def compute_ekin(veloc, masses):
    e = 0.0
    for i, vel in enumerate(veloc):
        e += np.sum(0.5*(vel**2)*masses[i])
    return e


def is_normal_mode_format(modes, natoms):
    """Check if the normal mode is in the given format!
        change the docstring to something useful!
    """
    import numpy as np

    thresh = 0.05
    nmodes = len(modes)
    #
    matrix = np.array([[mode.displacements[iatom][ixyz] for mode in modes]
                      for iatom in range(natoms)
                      for ixyz in range(3)])
    #
    result = np.dot(matrix.T, matrix)  # returns modes*modes matrix
    # compute the trace
    trace = np.trace(result)
    # set diagonal elemnts to 0.0
    np.fill_diagonal(result, 0.0)
    #
    nmodes_m1 = float(nmodes-1)
    # check that trace == nmodes - 1
    if trace > nmodes_m1:
        if abs((trace/nmodes) - 1.0) < thresh:
            # compute max/min ofdiagonal elements
            max_val = abs(np.max(result))
            min_val = abs(np.min(result))
            # check that there are no significant elements!
            if max_val < thresh and min_val < thresh:
                return True
    return False


def compute_norm_mode(mode, molecule):
    """Compute the norm of a mode"""
    import numpy as np
    norm = 0.0
    for iatom, displacement in enumerate(mode.displacements):
        for xyz in displacement:
            norm += xyz**2 * molecule.masses[iatom]/U_TO_AMU
    return np.sqrt(norm)


def scale_mode(imode, modes, expression):
    """example expression:

    lambda imode, iatom, ixyz, disp: disp/(norm/np.sqrt(molecule.masses[iatom]))

    """
    mode = deepcopy(modes[imode])
    #
    for iatom, displacement in enumerate(mode.displacements):
        mode.displacements[iatom] = [expression(imode, iatom, ixyz, disp)
                                     for ixyz, disp in enumerate(displacement)]
    #
    return mode


def create_mass_weighted_normal_modes(modes, molecule):
    """decide which format the normal modes are and convert to the mass-weighted once"""
    import numpy as np
    ANG_TO_BOHR = 1./0.529177211
    # compute norms
    norms = [compute_norm_mode(mode, molecule) for mode in modes]
    # define options
    options = {
            'gaussian-type': lambda imode, iatom, ixyz, disp: (disp/(norms[imode]
                                                               / np.sqrt(molecule.masses[iatom]
                                                               / U_TO_AMU))),
            'molpro-type': lambda imode, iatom, ixyz, disp: (disp * np.sqrt(
                                                             molecule.masses[iatom]
                                                             / U_TO_AMU)),
            'columbus-type': lambda imode, iatom, ixyz, disp: (disp * np.sqrt(
                                                               (molecule.masses[iatom]
                                                                / U_TO_AMU)
                                                               / ANG_TO_BOHR)),
            'mass-wighted': lambda imode, iatom, ixyz, disp: disp,
    }
    #
    possible_formats = []
    #
    for typ, expression in options.items():
        current_modes = [scale_mode(imode, modes, expression) for imode in range(len(modes))]
        if is_normal_mode_format(current_modes, molecule.natoms):
            selected_format = typ
            final_modes = current_modes
            possible_formats.append(selected_format)
    #
    if len(possible_formats) != 1:
        print("possible_formats = ", possible_formats)
        raise Exception('Could not specify format possible formats = "%s"'
                        % ", ".join(possible_formats))
    #
    return final_modes


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
        # paper paper says, that freq_factor is sqrt(2*PI*freq)
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


def wigner_gs(Q, P):
    """for a one-dimensional harmonic oscillator.
       Q contains the dimensionless coordinate of the
       oscillator and P contains the corresponding momentum."""
    return np.exp(-Q**2.0) * np.exp(-P**2.0)


def create_initial_conditions(settings, filename, molecule, modes, E_zero, model=False):
    # Initialize database
    if os.path.exists(filename):
        raise Exception('database "%s" does already exist' % filename)
    # create database
    db = Database(filename, settings)
    # set constants
    if model is False:
        db.set('atomids', molecule.atomids)
    db.set('masses', molecule.masses/U_TO_AMU)
    db.set('modes', np.array([mode.displacements for mode in modes]))
    db.set('freqs', np.array([mode.freq for mode in modes]))
    # add equilibrium values
    db.append('crd', molecule.crd)
    if model is False:
        db.append('veloc', np.zeros(molecule.natoms*3, dtype=np.double))
    else: 
        db.append('veloc',np.zeros(len(modes), dtype=np.double))
    db.append('energy', np.array([E_zero, E_zero, 0.0]))
    # go to next step
    db.increase
    #
    return db


if __name__ == '__main__':
    random.seed(16661)
    sampling = WignerSampling.from_molden('molden.dat')
    sampling.create_initial_conditions('init.db', 3, E_equil=0.0)
