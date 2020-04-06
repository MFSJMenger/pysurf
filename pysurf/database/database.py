from .dbtools import DatabaseRepresentation
from .dbtools import load_database


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
        nc = load_database(filename)
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
