from abc import abstractmethod
from collections import namedtuple
from copy import deepcopy
#
import numpy as np
#from pysurf.logger import get_logger
from ..utils.osutils import exists_and_isfile
from ..database import PySurfDB
from ..logger import Logger, get_logger
#
from ..system import Molecule
from ..spp import Model
from ..sampling.base_sampler import SamplerFactory
from .base_sampler import DynCondition, CrdCondition
from .normalmodes import Mode
#
from colt import Colt


class SamplingDB(PySurfDB):

    @classmethod
    def from_db(cls, dbfilename):
        """ Using an existing sampling file to set up sampling
            first information has to be read from database, then the database can be loaded.
            This is to make sure that databases don't get corrupted!
        """

        config = {'sampling_db': dbfilename}
        info = PySurfDB.info_database(dbfilename)
        config['n_conditions'] = info['length']
        #still a hard coded stuff as long as description of db not working
        config['method'] = 'Wigner'
        if 'atomids' in info['variables']:
            config['model'] = False
        else:
            config['model'] = True
        if 'veloc' in info['variables']:
            dynsampling = True
        else:
            dynsampling = False
        if info['dimensions']['frame'] == 1:
            sp = True
        else:
            sp = False
        db = cls.load_database(dbfilename, data=info['variables'], dimensions=info['dimensions'], model=config['model'], sp=sp)
        return db

    @classmethod
    def from_sampler(cls, config, sampler):
        getinit = sampler.get_init()

        system = getinit['system']
        modes = getinit['modes']
        if isinstance(system, Model):
            model = True
        else:
            model = False
        if sampler.condition.__name__ == "DynCondition":
            dynsampling = True
        else:
            dynsampling = False
        nmodes = len(modes)
        
        if model is False:
            natoms = system.natoms
            variables = ['model', 'crd_equi', 'masses', 'atomids', 'modes_equi', 'freqs_equi', 'crd', 'currstate']
            dimensions = {'natoms': natoms, 'nmodes': nmodes}
        else:
            variables = ['model', 'crd_equi', 'masses', 'modes_equi', 'freqs_equi', 'crd', 'currstate']
            dimensions = {'nmodes': nmodes}
        if dynsampling:
            variables += ['veloc']
        db = cls.generate_database(config['sampling_db'], data=variables, dimensions=dimensions, model=model)
        db.add_reference_entry(system, modes, model)
        return db

    @classmethod
    def create_db(cls, dbfilename, variables, dimensions, system, modes, model=False, sp=False):
        db = cls.generate_database(dbfilename, data=variables, dimensions=dimensions, model=model, sp=sp)
        db.add_reference_entry(system, modes, model)
        return db
    
    def append_condition(self, cond):
        self.append('crd', cond.crd)
        if self.dynsampling:
            self.append('veloc', cond.veloc)
            self.append('currstate', cond.state)
        self.increase

    def write_condition(self, cond, idx):
        self.set('crd', cond.crd, idx)
        if self.dynsampling:
            self.set('veloc', cond.veloc, idx)
            self.set('currstate', cond.state, idx)

    def write_xyz(self, filename, number):
        molecule = deepcopy(self.molecule)
        molecule.crd = self.get('crd', number)
        molecule.write_xyz(filename)

    def get_condition(self, idx):
        if idx >= self.nconditions:
            return None
        crd = self.get('crd', idx)
        if self.dynsampling:
            veloc = self.get('veloc', idx)
            state = int(self.get('currstate', idx))
            return self.condition(crd, veloc, state)
        return self.condition(crd)

    def get_config(self):
        config = {}
        config['n_conditions'] = self.nconditions
        config['method'] = 'Wigner'
        config['model'] = self.model
        return config

    @property
    def dynsampling(self):
        if 'veloc' in self.info['variables']:
            return True
        else:
            return False

    @property
    def condition(self):
        if self.dynsampling:
            return DynCondition
        else:
            return CrdCondition

    @property
    def nconditions(self):
        return len(self['crd'])

#    @classmethod
#    def _get_method_from_db(cls, db):
#        _method_number = db['method']
#        method = cls._get_method_from_number(_method_number)
#        return method
#
#    @staticmethod
#    def _get_method_from_number(numbers):
#        method = ''
#        for num in numbers:
#            if num != 0:
#                method += chr(num)
#            else:
#                break
#        return method
#
#    @staticmethod
#    def _get_number_from_method(method):
#        #length of numpy arry has to be consistent with length in settings
#        res = np.zeros(100, int)
#        for index, letter in enumerate(method):
#            res[index] = ord(letter)
#        return res
#
