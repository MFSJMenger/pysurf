import numpy as np

from pysurf.utils import exists_and_isfile
from .dbtools import DatabaseRepresentation, DatabaseGenerator
from .dbtools import load_database as l_db


class Database(object):
    """"Core Database, can store data given by the DatabaseRepresentation

    database is automatically build using a settings dictionary:

    dct = { 'dimensions': {
                 'frame': 'unlimited',
                 'natoms': 3,
                 'three': 3,
                 },
             'variables': {
                'crd': DBVariable(np.double, ('frame', 'natoms', 'three')),
             }
          }

    or string:

    dct = "

    [dims]
    frame = unlimited
    natoms = natoms
    three = 3
    [vars]
    crd = double :: (frame, natoms, three)

    """

    __slots__ = ('filename', '_rep', '_db', '_handle', '_closed', '_icurrent')

    def __init__(self, filename, settings, read_only=False):
        """Initialize new Database,
        if db exists:
           load existing database
           check that the settings of the old database
           are the same with the once used in the loading
           routine
        else:
           create new database
        """
        self.filename = filename
        #
        if isinstance(settings, dict):
            self._rep = DatabaseRepresentation(settings)
        elif isinstance(settings, str):
            self._rep = DatabaseRepresentation.from_string(settings)
        #
        self._closed = True
        #
        self._db, self._handle = self._rep.create_database(filename, False)
        #
        self._closed = False
        #
        self._icurrent = None

    @classmethod
    def load_db(cls, filename):
        nc = l_db(filename)
        rep = DatabaseRepresentation.from_db(nc)
        return cls(filename, {'variables': rep.variables, 'dimensions': rep.dimensions})

    @classmethod
    def empty_like(cls, filename, db):
        """ create an empty database with the same dimensions as db
            filename is the name of the new empty db """
        return cls(filename, {'variables': db._rep.variables, 'dimensions': db._rep.dimensions})

    def __getitem__(self, key):
        return self._handle.get(key, None)

    def __contains__(self, key):
        return key in self._handle

    def keys(self):
        return self._handle.keys()

    def get(self, key, ivalue):
        variable = self._handle[key]
        if variable.shape[0] > ivalue:
            return variable[ivalue]

    def get_keys(self):
        return self._db.variables.keys()

    def get_dimension_size(self, key):
        dim = self.db.dimensions.get(key, None)
        if dim is not None:
            return dim.size

    @property
    def closed(self):
        return self._closed

    @property
    def increase(self):
        self._icurrent += 1

    @property
    def info(self):
        return {'variables': list(self.get_keys()), 'dimensions': dict(self._rep.dimensions)}

    def append(self, key, value):
        """Append only for unlimited variables!"""
        variable = self._handle[key]
        unlimited = variable.get_dims()[0]
        assert(unlimited.isunlimited())
        if self._icurrent is None:
            self._icurrent = unlimited.size
        variable[self._icurrent, :] = value

    def set(self, key, value, ivalue=None):
        """set a given variable"""
        variable = self._handle[key]
        if variable.get_dims()[0].isunlimited():
            if ivalue is None:
                self.append(key, value)
            else:
                variable[ivalue, :] = value
        else:
            variable[:] = value

    def __del__(self):
        if self._closed is False:
            self._db.close()


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
        self.set('model', model, 0)
        if model is False:
            self.set('atomids', molecule.atomids, 0)
        self.set('masses', molecule.masses, 0)
        self.set('modes_equi', np.array([mode.displacements for mode in modes]), 0)
        self.set('freqs_equi', np.array([mode.freq for mode in modes]), 0)
        # add equilibrium values
        self.set('crd_equi', molecule.crd, 0)

    @property
    def masses(self):
        return np.array(self['masses'])
