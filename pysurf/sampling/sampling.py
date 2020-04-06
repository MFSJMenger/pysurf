from abc import abstractmethod
from copy import deepcopy
#
import numpy as np
#
from pysurf.logger import get_logger
from pysurf.utils.osutils import exists_and_isfile
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from .base_sampling import Condition
from pysurf.colt import PluginBase
#
from pysurf.molecule.molecule import Molecule
from .normalmodes import Mode


class SamplingFactory(PluginBase):

    _plugins_storage = '_methods'
    _is_plugin_factory = True

    _questions = """
    # chose the sampling method
    method =
    # database file where the conditions are saved or from which the conditions
    # are taken if the file already exists.
    filename = sampling.db
    # State whether the system is a model system
    model = False :: bool
    # Number of conditions that have to be created.
    number of conditions = 100 :: int
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                 for name, method in cls._methods.items()})

    @classmethod
    def from_inputfile(cls, inputfile):
        """ Class to create conditions due to user input. Initial conditions are saved
            in a file for further usage.
        """
        # Generate the config
        quests = cls.generate_questions('SAMPLING', config=None)
        config = quests.check_only(inputfile)
        return cls(config)

    @classmethod
    def from_db(cls, dbfilename):
        config = {'filename': dbfilename}
        db = Database.load_db(dbfilename)
        config['number of conditions'] = len(db['crd'])
        config['method'] = cls._get_method_from_db(db)
        if 'atomids' in db.keys():
            config['model'] = False
        else:
            config['model'] = True
        return cls(config)


class SamplingBase(SamplingFactory):
    """Basic sampling class"""

    subquestions: "inherited"
    #
    _register_plugin = False
    _questions = "inherited"
    #
    Condition = Condition

    def __init__(self, config, logger=None):
        """Setup sampling routine"""

        if logger is None:
            self.logger = get_logger('sampling.log', 'sampling')
        else:
            self.logger = logger
        # init data
        self._molecule = None
        self._modes = None
        self._equilibrium = None
        self._model = None
        # get sampling method
        self.method, self._method_number = self._get_method(config['method'][0])
        # setup sampler
        self.sampler = self._get_sampler(config['method'])
        # setup database
        self._setup_database(config)
        # check for conditions
        if self.nconditions < config['number of conditions']:
            self.logger.info('Adding additional entries to database of conditions')
            #self.add_conditions(config['number of conditions'] - self.nconditions)
            #self.nconditions = config['number of conditions']

    def _setup_database(self, config):
        """Load database, if it does not exist, create a new one"""
        # setup database
        if exists_and_isfile(config['filename']):
            self.logger.info(f"Reading conditions from file '{config['filename']}'")
            self._db = Database.load_db(config['filename'])
            self.nconditions = len(self._db['crd'])
        else:
            self.logger.info(f"Setting up new database for conditions: {config['filename']}")
            self.nconditions = 0
            self._model = config['model']
            self._setup_db(config)
        # sanity check
        if self.nconditions > config['number of conditions']:
            self.logger.error("Number of requested conditions less than current conditions!")

    def export_condition(self, dbfilename, number):
        """Export conditions"""
        db = Database.empty_like(dbfilename, self._db)
        self._add_reference_entry(db, self.molecule, self.modes, self.method, self.model)
        self._write_condition(db, self.get_condition(number))

    def write_xyz(self, filename, number):
        molecule = deepcopy(self.molecule)
        molecule.crd = self._db.get('crd', number)
        molecule.write_xyz(filename)

    def __iter__(self):
        self._start = 1  # skip equilibrium structure
        return self

    def __next__(self):
        cond = self.get_condition(self._start)
        if cond is not None:
            self._start += 1
            return cond
        raise StopIteration

    def get_condition(self, idx):
        if idx >= self.nconditions:
            return None
        crd = self._db.get('crd', idx)
        return self.Condition(crd)

    def add_conditions(self, nconditions, state=0):
        # TODO
        # Take the random seed and random counter from the database to
        # assure the consistency with all
        # the previous conditions
        for _ in range(nconditions):
            cond = self.sampler.get_condition()
            self._write_condition(self._db, cond)

    def _setup_db(self, config):
        getinit = self.sampler.get_init()
        molecule = getinit['molecule']
        self._modes = getinit['modes']
        self.natoms = molecule.natoms
        self.nmodes = len(self.modes)
        self._db = Database(config['filename'], self._settings)
        self._add_reference_entry(self._db, molecule, self.modes, self.method, self.model)

    @property
    def equilibrium(self):
        if self._equilibrium is None:
            self._equilibrium = self.get_condition(0)
        return self._equilibrium

    @equilibrium.setter
    def equilibrium(self, value):
        self._equilibrium = value

    @property
    def molecule(self):
        if self._molecule is None:
            self._molecule = Molecule(np.copy(self._db['atomids']),
                                      np.copy(self._db.get('crd', 0)),
                                      np.copy(self._db['masses']))
        return self._molecule

    @molecule.setter
    def molecule(self, value):
        self._molecule = value

    @property
    def modes(self):
        if self._modes is None:
            self._modes = [Mode(freq, mode) for freq, mode in zip(self._db['freqs'],
                                                                  self._db['modes'])]
        return self._modes

    @modes.setter
    def modes(self, value):
        self._modes = value

    @property
    def model(self):
        if self._model is None:
            if 'atomids' in self._db.keys():
                self._model = False
            else:
                self._model = True
        return self._model

    @model.setter
    def model(self, value):
        self._model = value

    @property
    def _settings(self):
        if self.model is False:
            settings = {
                    'dimensions': {
                        'frame': 'unlimited',
                        'nmodes': self.nmodes,
                        'natoms': self.natoms,
                        'three': 3,
                        'one': 1,
                        },
                    'variables': {
                        'modes':  DBVariable(np.double, ('nmodes', 'natoms', 'three')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'atomids': DBVariable(np.integer, ('natoms',)),
                        'masses': DBVariable(np.double, ('natoms',)),
                        'method': DBVariable(np.integer, ('one',)),
                        'crd': DBVariable(np.double, ('frame', 'natoms', 'three')),
                    },
            }

        if self.model is True:
            settings = {
                    'dimensions': {
                        'frame': 'unlimited',
                        'nmodes': self.nmodes,
                        'three': 3,
                        'one': 1,
                        },
                    'variables': {
                        'modes':  DBVariable(np.double, ('nmodes', 'nmodes')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'masses': DBVariable(np.double, ('nmodes',)),
                        'method': DBVariable(np.integer, ('one',)),
                        'crd': DBVariable(np.double, ('frame', 'nmodes')),
                    },
            }
        return settings

    def _get_sampler(self, config):
        return self._methods[self.method].from_config(config)

    def _get_method(self, method):
        number = tuple(self._methods.keys()).index(method)
        return method, number

    @classmethod
    def _get_method_from_db(cls, db):
        _method_number = db['method'][0]
        method = cls._get_method_from_number(_method_number)
        return method, _method_number

    @staticmethod
    def _get_method_from_number(number):
        try:
            out = tuple(SamplingBase._methods.keys())[number]
        except IndexError:
            raise Exception("Sampling method not known!") from None
        return out

    @staticmethod
    def _get_number_from_method(method):
        return list(SamplingBase._methods.keys()).index(method)

    @staticmethod
    def _add_reference_entry(db, molecule, modes, method, model):
        if model is False:
            db.set('atomids', molecule.atomids)
        db.set('masses', molecule.masses)
        db.set('modes', np.array([mode.displacements for mode in modes]))
        db.set('freqs', np.array([mode.freq for mode in modes]))
        _method_number = SamplingBase._get_number_from_method(method)
        db.set('method', _method_number)
        # add equilibrium values
        db.append('crd', molecule.crd)
        # go to next step
        db.increase

    @staticmethod
    def _write_condition(db, cond):
        db.append('crd', cond.crd)
        db.increase
