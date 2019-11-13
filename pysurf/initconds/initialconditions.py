import os
from copy import deepcopy
from collections import namedtuple, OrderedDict
#
import numpy as np
import numpy.random as random
#
from pysurf.molden import MoldenParser
#
from pysurf.logger import get_logger
from pysurf.utils.osutils import exists_and_isfile
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import DatabaseTools
from wigner import WignerSampling
#
from pysurf.molecule.molecule import Molecule
#
from pysurf.colt import AskQuestions
from pysurf.colt import QuestionGenerator
#  used to save just InitialConditions (coordinates, velocities) for a given molecule


class InitialConditions(object):
    _methods = OrderedDict({'wigner': WignerSampling})
    questions = """
    # database file where the initial conditions are saved or from which the initial conditions
    # are taken if the file already exists.
    outputfile = initconds.db


    # Describes which sampling algorithm is used to generate the initial conditions.
    # The default is wigner.
    sampling = wigner

    # State whether the system you are describing is a model system
    model = False :: bool

    # Initial state of the trajectories
    initial state = 0 :: int

    # Number of initial conditions that have to be created.
    # The default value is 100.
    number of initial conditions = 100 :: int
    """

    def __init__(self, inputfile, logger=None):
        """ Class to create initial conditions due to user input. Initial conditions are saved 
            in a file for further usage.
        """
        # Generate the config
        quests = QuestionGenerator(self.questions)
        quests.generate_cases("", "sampling", {
            name: sampling.questions for name, sampling in self._methods.items()})
        quests = AskQuestions("INITIAL CONDITIONS", quests.questions)
        self.config = quests.ask(inputfile)
        

        if logger is None:
            self.logger = get_logger('initconds.log', 'initconds')
        else:
            self.logger = logger

        

        if exists_and_isfile(self.config['outputfile']):
            self._db = Database.load_db(self.config['outputfile'])
            self.nconditions = len(self._db['crd'])
            self._get_method_from_db()
        else:
            self.nconditions = 0
            self.model = self.config['model']
            self._setup_method()
            sampler = self._methods[self.method](self.config)
            self._setup_db()
            self._add_reference_entry()


#    @classmethod
#    def from_conditions(cls, filename, molecule, modes, nconditions=1000, E_equil=0.0, model=False):
#        db = cls.create_initial_conditions(filename, molecule, modes, nconditions, E_equil, model)
#        return cls(db, nconditions)
#
#    @classmethod
#    def save_condition(cls, filename, init_cond):
#        natoms = len(init_cond.crd)
#        nmodes = natoms*3-6
#        db = Database(filename, cls.generate_settings(natoms=natoms, nmodes=nmodes, model=False))
#        db.append('veloc', init_cond.veloc)
#        db.append('crd', init_cond.crd)
#
#    @classmethod
#    def from_db(cls, filename, E_quil=0.0):
#        db_settings = DatabaseTools.get_variables(filename, ["natoms", "nmodes"])
#        db = Database(filename, cls.generate_settings(db_settings['natoms'], 
#                                                      db_settings['nmodes']))
#        return cls(db, db['crd'].shape[0])

    def _setup_db(self):
        if self.model is False:
            settings = {
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
                        'state': DBVariable(np.double, ('frame', 'one'))
                    },
            }

        if self.model is True:
            settings = {
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
                        'state': DBVariable(int, ('frame', 'one'))
                    },
            }

        filename = self.config['outputfile']
        self._db = Database(settings, filename)

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

        sampler = self._methods[self.method](self.config)

        for _ in range(nconditions):
            cond = sampler.get_condition()
            self._write_initial_condition(cond)



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

#    @classmethod
#    def create_initial_conditions(cls, filename, molecule, modes, nconditions, E_equil, state=0, model=False):
#        """ """
#        if model is False:
#            db = create_initial_conditions(cls.generate_settings(molecule.natoms, len(modes)),
#                                       filename, molecule, modes, E_equil, model)
#        else:
#            db = create_initial_conditions(cls.generate_settings(nmodes=len(modes), model=True),
#                                       filename, molecule, modes, E_equil, model)
#
#        for _ in range(nconditions):
#            e_pot, conds = get_initial_condition(molecule, modes)
#            db.append('crd', conds.crd)
#            db.append('veloc', conds.veloc)
#            db.append('state', state)
#            e_kin = compute_ekin(conds.veloc, molecule.masses)
#            db.append('energy', np.array([e_kin+e_pot, e_pot, e_kin]))
#            db.increase
#        return db

    def _get_method_from_db():
        self._method_number = self._db['method']
        self.method = self._get_method_from_number()

    def _get_method_from_number(self, number):
        return tuple(self._methods.keys())[number]

    def _get_number_from_method(self, method):
        return list(self._methods.keys()).index(self.method)

    def _setup_method(self):
        self.method = self.config['sampling']
        self._method_number = self._get_number_from_method(self.method)
        
    def _add_reference_entry():

        init = self.sampler.get_init()
        molecule = init['molecule']
        modes = init['modes']

        if self.model is False:
            self._db.set('atomids', molecule.atomids)
        self._db.set('masses', molecule.masses/U_TO_AMU)
        self._db.set('modes', np.array([mode.displacements for mode in modes]))
        self._db.set('freqs', np.array([mode.freq for mode in modes]))
        self._db.set('method', self._method_number)
        # add equilibrium values
        self._db.append('crd', molecule.crd)
        if self.model is False:
            self._db.append('veloc', np.zeros(molecule.natoms*3, dtype=np.double))
        else: 
            self._db.append('veloc',np.zeros(len(modes), dtype=np.double))
        self._db.append('energy', np.array([E_zero, E_zero, 0.0]))
        # go to next step
        self._db.increase

    def _write_initial_condition(self, cond):
        self._db.append('crd', cond.crd)
        self._db.append('veloc', cond.veloc)
        self._db.append('state', state)
        self._db.increase



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
