import os
from copy import deepcopy
from collections import namedtuple
#
import numpy as np
import numpy.random as random
#
from pysurf.molden import MoldenParser
#
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import DatabaseTools
#
from pysurf.molecule.molecule import Molecule
#
#  used to save just InitialConditions (coordinates, velocities) for a given molecule
InitialCondition = namedtuple("InitialCondition", ['crd', 'veloc'])


class InitialConditions(object):
    def __init__(self, inputfile, logger=None):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """

        question_string = """
        # database file where the initial conditions are saved or from which the initial conditions
        # are taken if the file already exists.
        outputfile = initconds.db

        # Input source for the normal modes and/or frequencies, which are used to generate the 
        # initial conditions.
        # Possible options are:
        # - molden
        # - frequencies
        from = none :: str :: [molden, frequencies]

        # Describes which sampling algorithm is used to generate the initial conditions.
        # The default is wigner.
        sampling = wigner

        # Number of initial conditions that have to be created.
        # The default value is 100.
        number of initial conditions = 100 :: int

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

        self.logger = get_logger('initconds.log', 'initconds')

        quests = AskQuestions.from_string("INITIAL CONDITIONS", question_string, config=inputfile)
        if logger is None:
            self.logger = get_logger('initconds.log', 'initconds')
        else:
            self.logger = logger

        if exists_and_isfile(filename):
            self._db = load_db(filename)
        self.nconditions = nconditions
        self._db = db
        self._molecule = None

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

    def add_initial_conditions(self, nconditions, state=0):
        # TODO 
        # Take the random seed and random counter from the database to assure the consistency with all
        # the previous conditions

        modes = [Mode(freq, mode) for freq, mode in zip(self._db['freqs'], self._db['modes'])]

        for _ in range(nconditions):
            e_pot, conds = get_initial_condition(self.molecule, modes)
            self._db.append('crd', conds.crd)
            self._db.append('veloc', conds.veloc)
            self._db.append('state', state)
            e_kin = compute_ekin(conds.veloc, self.molecule.masses)
            self._db.append('energy', np.array([e_kin+e_pot, e_pot, e_kin]))
            self._db.increase


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
        if self._molecule is None:
            self._molecule = Molecule(np.copy(self._db['atomids']),
                                      np.copy(self._db.get('crd', 0)),
                                      np.copy(self._db['masses']))
        return self._molecule

    @classmethod
    def create_initial_conditions(cls, filename, molecule, modes, nconditions, E_equil, state=0, model=False):
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
            db.append('state', state)
            e_kin = compute_ekin(conds.veloc, molecule.masses)
            db.append('energy', np.array([e_kin+e_pot, e_pot, e_kin]))
            db.increase
        return db





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
