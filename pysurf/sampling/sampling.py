from abc import abstractmethod
from collections import namedtuple
from copy import deepcopy
#
import numpy as np
#from pysurf.logger import get_logger
from pysurf.utils.osutils import exists_and_isfile
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.colt import Colt
from pysurf.logger import Logger, get_logger
#
from pysurf.molecule.molecule import Molecule
from pysurf.sampling.base_sampler import SamplerFactory
from .base_sampler import DynCondition, CrdCondition
from .normalmodes import Mode




class Sampling(Colt):
    """ Sampling is the header class for the sampling routines. It asks the main questions, selects
        the sampler and reads and writes the conditions to the sampling database.
    """
    _questions = """
    # Number of Calculations
    n_conditions = 100 :: int

    method = :: str

    # Database containing all the initial conditions
    sampling_db = sampling.db :: file

    #State whether the system is a model system
    model = False :: bool
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("method", {name: method.questions
                                 for name, method in SamplerFactory._methods.items()})

    @classmethod
    def from_config(cls, config):
        if exists_and_isfile(config['sampling_db']):
            db = SamplingDB.from_db(config['sampling_db'])
            sampler = None
        else:
            sampler = cls._get_sampler(config['method'])
            db = SamplingDB.create_db(config, sampler)
        return cls(config, db, sampler)

    @classmethod
    def from_inputfile(cls, inputfile):
        # Generate the config
        quests = cls.generate_questions('SAMPLING', config=None)
        config = quests.check_only(inputfile)
        return cls.from_config(config)    


    @classmethod
    def from_db(cls, dbfilename):
        if exists_and_isfile(dbfilename):
            db = SamplingDB.from_db(dbfilename)
            config = db.get_config()
            return cls(config, db=db)
        else:
            logger = get_logger('sampling.log', 'sampling')
            logger.error(f'File of sampling database {dbfilename} not found')

    def __init__(self, config, db, sampler=None, logger=None):
        """ Sampling always goes with a database, if not needed use Sampler class
        """

        self._db = db
        self.config = config
        self.sampler = sampler

        if logger is None:
            self.logger = get_logger('sampling.log', 'sampling')
        else:
            self.logger = logger

        # check for conditions
        if self._db.nconditions < config['n_conditions']:
            # setup sampler
            if sampler is None:
                self.sampler = self._get_sampler(config['method'])
            else:
                self.sampler = sampler
            self.logger.info('Adding additional entries to database of conditions')
            self.add_conditions(config['n_conditions'] - self.nconditions)

    def add_conditions(self, nconditions, state=0):
        # TODO
        # Take the random seed and random counter from the database to
        # assure the consistency with all
        # the previous conditions
        for _ in range(nconditions):
            cond = self.sampler.get_condition()
            self._db.write_condition(cond)

    def get_condition(self, idx):
        return self._db.get_condition(idx)

    def write_xyz(self, filename, idx):
        return self._db.write_xyz(filename, idx)

    def __iter__(self):
        self._start = 1  # skip equilibrium structure
        return self

    def __next__(self):
        cond = self.get_condition(self._start)
        if cond is not None:
            self._start += 1
            return cond
        raise StopIteration

    def export_condition(self, dbnew, number):
        return self._db.export_condition(dbnew, number)

    @staticmethod
    def _get_sampler(config):
        return SamplerFactory._methods[config.value].from_config(config)
    
    @property
    def nconditions(self):
        return self._db.nconditions

    @property
    def equilibrium(self):
        return self._db.equilibrium

    

class SamplingDB:

    @classmethod
    def from_db(cls, dbfilename):
        config = {'sampling_db': dbfilename}
        db = Database.load_db(dbfilename)
        config['n_conditions'] = len(db['crd'])
        config['method'] = cls._get_method_from_db(db)
        if 'atomids' in db.keys():
            config['model'] = False
        else:
            config['model'] = True
        if 'veloc' in db.keys():
            dynsampling = True
        else:
            dynsampling = False
        return cls(config, db, dynsampling)

    @classmethod
    def create_db(cls, config, sampler):
        getinit = sampler.get_init()
        molecule = getinit['molecule']
        modes = getinit['modes']
        model = config['model']
        if sampler.condition.__name__ == "DynCondition":
            dynsampling = True
        else:
            dynsampling = False
        natoms = molecule.natoms
        nmodes = len(modes)
        settings = cls._get_settings(natoms=natoms, nmodes=nmodes, dynsampling=dynsampling, model=model)
        db = Database(config['sampling_db'], settings)
        cls._add_reference_entry(db, molecule, modes, config['method'].value, model)
        return cls(config, db, dynsampling)

    def __init__(self, config, db, dynsampling):
        self._db = db
        self.config = config
        self._molecule = None
        self._modes = None
        self._model = None
        self.dynsampling = dynsampling
        if self.dynsampling:
            self.condition = DynCondition
        else:
            self.condition = CrdCondition

    def export_condition(self, dbfilename, number):
        """Export conditions"""
        db = Database.empty_like(dbfilename, self._db)
        self._add_reference_entry(db, self.molecule, self.modes, self.method, self.model)
        samplingdb = SamplingDB(self.config, db, self.dynsampling)
        samplingdb.write_condition(self.get_condition(number))
        return samplingdb

    def get_config(self):
        return self.config

    def write_condition(self, cond):
        self._db.append('crd', cond.crd)
        if self.dynsampling:
            self._db.append('veloc', cond.veloc)
            self._db.append('state', cond.state)
        self._db.increase

    def write_xyz(self, filename, number):
        molecule = deepcopy(self.molecule)
        molecule.crd = self._db.get('crd', number)
        molecule.write_xyz(filename)

    def get_condition(self, idx):
        if idx >= self.nconditions:
            return None
        crd = self._db.get('crd', idx)
        if self.dynsampling:
            veloc = self._db.get('veloc', idx)
            #state = self._db.get('state', idx)
            return self.condition(crd, veloc, None)
        return self.condition(crd)

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
    def method(self):
        if isinstance(self.config['method'], str):
            return self.config['method']
        else:
            return self.config['method'].value
            
    @method.setter
    def method(self, value):
        self.method = value

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
    def nconditions(self):
        return len(self._db['crd'])


    @staticmethod
    def _get_settings(nmodes, natoms, dynsampling, model):
        if model is False:
            settings = {
                    'dimensions': {
                        'frame': 'unlimited',
                        'nmodes': nmodes,
                        'natoms': natoms,
                        'three': 3,
                        'one': 1,
                        'hundred': 100,
                        },
                    'variables': {
                        'modes':  DBVariable(np.double, ('nmodes', 'natoms', 'three')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'atomids': DBVariable(np.integer, ('natoms',)),
                        'masses': DBVariable(np.double, ('natoms',)),
                        'method': DBVariable(np.integer, ('hundred',)),
                        'crd': DBVariable(np.double, ('frame', 'natoms', 'three')),
                    },
            }
            if dynsampling:
                settings['variables']['veloc'] = DBVariable(np.double, ('frame', 'natoms', 'three'))
        if model is True:
            settings = {
                    'dimensions': {
                        'frame': 'unlimited',
                        'nmodes': self.nmodes,
                        'three': 3,
                        'one': 1,
                        'hundred': 100,
                        },
                    'variables': {
                        'modes':  DBVariable(np.double, ('nmodes', 'nmodes')),
                        'freqs':  DBVariable(np.double, ('nmodes',)),
                        'masses': DBVariable(np.double, ('nmodes',)),
                        'method': DBVariable(np.integer, ('hundred',)),
                        'crd': DBVariable(np.double, ('frame', 'nmodes')),
                    },
            }

            if dynsampling:
                settings['variables']['veloc'] = DBVariable(np.double, ('frame', 'nmodes'))
        if dynsampling:
            settings['variables']['state'] = DBVariable(np.double, ('frame', 'one'))
        return settings

    @classmethod
    def _get_method_from_db(cls, db):
        _method_number = db['method']
        method = cls._get_method_from_number(_method_number)
        return method

    @staticmethod
    def _get_method_from_number(numbers):
        method = ''
        for num in numbers:
            if num != 0:
                method += chr(num)
            else:
                break
        return method

    @staticmethod
    def _get_number_from_method(method):
        #length of numpy arry has to be consistent with length in settings
        res = np.zeros(100, int)
        for index, letter in enumerate(method):
            res[index] = ord(letter)
        return res

    @staticmethod
    def _add_reference_entry(db, molecule, modes, method, model):
        if model is False:
            db.set('atomids', molecule.atomids)
        db.set('masses', molecule.masses)
        db.set('modes', np.array([mode.displacements for mode in modes]))
        db.set('freqs', np.array([mode.freq for mode in modes]))
        _method_number = SamplingDB._get_number_from_method(method)
        db.set('method', _method_number)
        # add equilibrium values
        db.append('crd', molecule.crd)
        # go to next step
        db.increase


