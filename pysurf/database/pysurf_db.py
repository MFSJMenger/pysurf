from collections import namedtuple
import numpy as np

from pysurf.utils import exists_and_isfile
from .dbtools import DatabaseRepresentation, DatabaseGenerator
from .dbtools import load_database as l_db
from .database import Database
from pysurf.molecule import Molecule, Mode



class PySurfDB(Database):

    _dimensions = {
            'frame': 'unlimited',
            'natoms': None,
            'nstates': None,
            'nmodes': None,
            'three': 3,
            'two': 2,
            'one': 1,
    }
    _variables = DatabaseGenerator("""
        [variables]
        crd_equi  = double :: (natoms, three)
        atomids   = int    :: (natoms)
        freqs_equi= double :: (nmodes)
        modes_equi= double :: (nmodes, natoms, three)
        masses    = double :: (natoms)
        model     = int    :: (one)

        crd       = double :: (frame, natoms, three)
        veloc     = double :: (frame, natoms, three)
        accel     = double :: (frame, natoms, three)
        energy    = double :: (frame, nstates)
        gradient  = double :: (frame, nstates, natoms, three)
        fosc      = double :: (frame, nstates)
        transmom  = double :: (frame, nstates, three)
        currstate = double :: (frame, one)
        ekin      = double :: (frame, one)
        epot      = double :: (frame, one)
        etot      = double :: (frame, one)
        nacs      = double :: (frame, nstates, nstates, natoms, three)
    """)['variables']


    @classmethod
    def _get_settings(cls, variables):
        out = {'variables': {}, 'dimensions': {}}
        varis = out['variables']
        dims = out['dimensions']
       
        for var in variables:
            try:
                varis[var] = cls._variables[var]
            except:
                raise Exception(f"Variable {var} unknown")
            for dim in varis[var].dimensions:
                try:
                    dims[dim] = cls._dimensions[dim]
                except:
                    raise Exception(f"Dimension {dim} unknown")
        return out

    @classmethod
    def _prepare_settings(cls, settings, dimensions, model, single_point):
        dims = settings['dimensions']
        for dim, value in settings['dimensions'].items():
            if value is None:
                dims[dim] = dimensions[dim]
            if value == 'unlimited' and single_point is True:
                dims[dim] = 1
        # TODO: modify dimensions for the case db, etc
    
    @classmethod
    def info_database(cls, filename):
        db = Database.load_db(filename)
        info = {'variables':[]}
        for var in cls._variables.keys():
            if var in db:
                info['variables'] += [var]
        info['dimensions'] = db._rep.dimensions
        info['length'] = len(db['crd'])
        return info

    @classmethod
    def generate_database(cls, filename, data=None, dimensions=None, units=None, attributes=None, descriptition=None, model=False, sp=False):
        if dimensions is None:
            dimensions = {}
        if data is None:
            data = []
        settings = cls._get_settings(data)
        cls._prepare_settings(settings, dimensions, model, sp)
        return cls(filename, settings)

    @classmethod
    def load_database(cls, filename, data=None, dimensions=None, units=None, attributes=None, descriptition=None, model=False, sp=False, read_only=False):
        if read_only is True:
            return Database.load_db(filename)
        #
        if not exists_and_isfile(filename):
            raise Exception(f"Cannot load database {filename}")
        return cls.generate_database(filename, data, dimensions, units, attributes, descriptition, model, sp)

    def add_reference_entry(self, molecule, modes, model):
        self.set('model', int(model))
        if model is False:
            self.set('atomids', molecule.atomids)
        self.set('masses', molecule.masses)
        self.set('modes_equi', np.array([mode.displacements for mode in modes]))
        self.set('freqs_equi', np.array([mode.freq for mode in modes]))
        # add equilibrium values
        self.set('crd_equi', molecule.crd)

    @property
    def masses(self):
        return np.array(self['masses'])

    @property
    def modes(self):
        self._modes = [Mode(freq, mode) for freq, mode in zip(np.copy(self['freqs_equi']), np.copy(self['modes_equi']))]
        return self._modes


    @property
    def molecule(self):
        self._molecule = Molecule(np.copy(self['atomids']),
                                      np.copy(self['crd_equi']),
                                      np.copy(self['masses']))
        return self._molecule 

    @property
    def model(self):
            return bool(self['model'])
