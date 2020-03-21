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
from .base_sampling import Condition
from pysurf.colt import PluginBase
#
from pysurf.molecule.molecule import Molecule
from .normalmodes import Mode 
#
from pysurf.colt import Colt

class SamplingFactory(PluginBase):
    _plugins_storage = '_methods'
    _is_plugin_factory = True

    _questions = """
    method = 
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
        quests = cls.generate_questions('SAMPLING', config=inputfile)
        config = quests.ask(inputfile)
        return cls(config)


    @classmethod
    def from_db(cls, dbfilename):
        config = {'outputfile': dbfilename}
        db = Database.load_db(dbfilename)
        config['number of conditions'] = len(db['crd'])
        config['method'] = SamplingBase._get_method_from_db(db)
        if 'atomids' in db.keys():
            config['model'] = False
        else:
            config['model'] = True
        return cls(config)


class SamplingBase(SamplingFactory):
    _register_plugin = False
    subquestions : 'inherited'
    _questions = 'inherited'

    def _generate_subquestions(cls, questions):
        questions.add_questions_to_block("""
             # database file where the conditions are saved or from which the conditions
             # are taken if the file already exists.
             outputfile = sampling.db

             # State whether the system is a model system
             model = False :: bool

             # State whether the system you are describing is a model system
             # Number of conditions that have to be created.
             number of conditions = 100 :: int
             """)


    def __init__(self, config, logger = None):
        self.config = config

        if logger is None:
            self.logger = get_logger('sampling.log', 'sampling')
        else:
            self.logger = logger

        self._molecule = None 
        self._modes = None
        self._equilibrium = None
        self._model = None
        self.condition = Condition

        if exists_and_isfile(self.config['outputfile']):
            self.logger.info('Database of conditions exists')
            self.logger.info('Taking info from database ' + self.config['outputfile'])
            self._db = Database.load_db(self.config['outputfile'])
            self.nconditions = len(self._db['crd'])
            self.method, self._method_number = SamplingBase._get_method_from_db(self._db)
            if self.nconditions < self.config['number of conditions']:
                self.logger.info('Adding additional entries to database of conditions')
                self.add_conditions(self.config['number of conditions']
                                            - self.nconditions)
        else:
            self.logger.info('Setting up new database for conditions: ' 
                             + self.config['outputfile'])
            self.nconditions = 0
            self._model = self.config['model']
            self._setup_method()
            self._setup_db()
            self.logger.info('Adding {0} new conditions to database'.format(
                             self.config['number of conditions']))
            self.add_conditions(self.config['number of conditions'])
            self.nconditions = self.config['number of conditions']


    def export_condition(self, dbfilename, number):
        db = Database.empty_like(dbfilename, self._db)
        SamplingBase._add_reference_entry(db, self.molecule, self.modes, self.method, self.model)
        SamplingBase._write_condition(db, self.get_condition(number))        
    
        

    def __iter__(self):
        self._start = 1 # skip equilibrium structure
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
        return Condition(crd)

    def add_conditions(self, nconditions, state=0):
        # TODO 
        # Take the random seed and random counter from the database to assure the consistency with all
        # the previous conditions

        sampler = self._get_sampler()
        for _ in range(nconditions):
            cond = sampler.get_condition()
            self._write_condition(self._db, cond)

    def _setup_db(self):
        sampler = self._methods[self.method].from_config(self.config['method'])
        getinit = sampler.get_init()
        molecule = getinit['molecule']
        self._modes = getinit['modes']
        self.natoms = molecule.natoms
        self.nmodes = len(self.modes)
        self._db = Database(self.config['outputfile'], self._settings)
        SamplingBase._add_reference_entry(self._db, molecule, self.modes, self.method, self.model)


    @property
    def equilibrium(self):
        if self._equilibrium is None:
            return self.get_condition(0)
        else:
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
                        'modes':  DBVariable(np.double, ('nmodes','nmodes')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'masses': DBVariable(np.double, ('nmodes',)),
                        'method': DBVariable(np.integer, ('one',)),
                        'crd': DBVariable(np.double, ('frame', 'nmodes')),
                    },
            }
        return settings

    def _get_sampler(self):
        return self._methods[self.method].from_config(self.config['method'])

    @staticmethod
    def _get_method_from_db(db):
        _method_number = db['method'][0]
        method = SamplingBase._get_method_from_number(_method_number)
        return method, _method_number

    @staticmethod
    def _get_method_from_number(number):
        return tuple(SamplingBase._methods.keys())[number]
    
    @staticmethod
    def _get_number_from_method(method):
        return list(SamplingBase._methods.keys()).index(method)

    def _setup_method(self):
        self.method = self.config['method'].value
        self._method_number = SamplingBase._get_number_from_method(self.method)
        
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




