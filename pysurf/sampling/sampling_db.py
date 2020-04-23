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


class SamplingDB(Colt):

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
        cls._add_reference_entry(db, molecule, modes, config['method'].value, model, dynsampling)
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
        self._add_reference_entry(db, self.molecule, self.modes, self.method, self.model, self.dynsampling)
        config = self.config
        config['n_conditions'] = 1
        samplingdb = SamplingDB(config, db, self.dynsampling)
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

    def set_state(self, state, idx):
        if self.dynsampling:
            self._db.set('state', state, idx)

    def set_crd(self, crd, idx):
        self._db.set('crd', crd, idx)

    def set_veloc(self, veloc, idx):
        if self.dynsampling:
            self._db_set('veloc', veloc, idx)

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
            state = self._db.get('state', idx)
            return self.condition(crd, veloc, state)
        return self.condition(crd)

    def get_atomids(self):
        return self._db['atomids']

    def get_masses(self):
        return self._db['masses']

    @property
    def equilibrium(self):
        if self._equilibrium is None:
            self._equilibrium = self.get_condition(0)
        return self._equilibrium

    @equilibrium.setter
    def equilibrium(self, value):
        elf._equilibrium = value

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
    def natoms(self):
        return self.molecule.natoms

    @property
    def atomids(self):
        return self.molecule.atomids

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

    @property
    def mass(self):
        return self.molecule.masses

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
            settings['variables']['state'] = DBVariable(np.integer, ('frame', 'one'))
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
    def _add_reference_entry(db, molecule, modes, method, model, dynsampling):
        if model is False:
            db.set('atomids', molecule.atomids)
        db.set('masses', molecule.masses)
        db.set('modes', np.array([mode.displacements for mode in modes]))
        db.set('freqs', np.array([mode.freq for mode in modes]))
        _method_number = SamplingDB._get_number_from_method(method)
        db.set('method', _method_number)
        # add equilibrium values
        db.append('crd', molecule.crd)
        if dynsampling:
            db.append('veloc', np.zeros(molecule.crd.shape))
        # go to next step
        db.increase


