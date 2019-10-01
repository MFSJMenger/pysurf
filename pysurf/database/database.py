import netCDF4

from .dbtools import DatabaseRepresentation
from ..utils.osutils import exists_and_isfile


class Database(object):

    __slots__ = ('filename', '_rep', '_db', '_handle', 'closed')

    def __init__(self, filename, settings):
        """Initialize new Database,
        if db exisits:
           load existing database
           check that the settings of the old database
           are the same with the once used in the loading
           routine
        else:
           create new database
        """
        #
        self.filename = filename
        self._rep = DatabaseRepresentation(settings)
        # 
        if exists_and_isfile(filename):
            self._db, self._handle = self._load_database(filename)
        else:
            self._db, self._handle = self._init_database(filename)
        self.closed = False

    def __getattr__(self, key):
        if key in self._handle:
            return self._handle[key]
        return object.__getattribute__(self, key)

    def get_dimension_size(self, key):
        if key in self._db.dimensions.keys():
            return self._db.dimensions[key].size

    def append(self, key, value):
        variable = self._handle[key]
        unlimited = variable.get_dims()[0]
        assert(unlimited.isunlimited())
        variable[unlimited.size, :] = value

    def close(self):
        if self.closed is False:
            self._db.close()
            self.closed = True

    def __del__(self):
        self.close()

    def _init_database(self, filename):
        """Create a new database"""
        return create_dataset(filename, self._rep)

    def _load_database(self, filename):
        """Load an existing database and check
           that it is compatable with the existing one"""
        db = load_database(filename)
        ref = DatabaseRepresentation.from_db(db)
        if ref != self._rep:
            raise Exception('Database is not in agreement with ask settings!')
        return db, db.variables


def create_dataset(filename, settings):
    # 
    nc = netCDF4.Dataset(filename, 'w')
    # create dimensions
    for dim_name, dim in settings['dimensions'].items():
        if dim == 'unlimited':
            dim = None
        nc.createDimension(dim_name, dim)
    # create variables
    handle = {}
    for var_name, variable in settings['variables'].items():
        handle[var_name] = nc.createVariable(var_name, variable.type, variable.dimensions)
    return nc, handle 


def load_database(filename, io_options='a'):
    return netCDF4.Dataset(filename, io_options)
