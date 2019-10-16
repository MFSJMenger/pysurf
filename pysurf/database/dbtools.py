"""Tools to store information on the Variables and Dimensions in the Database"""
import netCDF4
import numpy as np

from ..utils.osutils import exists_and_isfile


class DBVariable(object):
    """Store info for Database and easy comparison"""

    __slots__ = ('type', 'dimensions')

    def __init__(self, typ, dim):
        self.type = typ
        self.dimensions = dim

    def __eq__(self, rhs):
        assert(isinstance(rhs, self.__class__))
        if len(self.dimensions) != len(rhs.dimensions):
            return False
        if self.type != rhs.type or any(self.dimensions[i] != rhs.dimensions[i] for i in range(len(self.dimensions))):
            print(self, rhs)
            return False
        return True

    def __str__(self):
        return f"DBVariable(type = {self.type}, dimension = {self.dimensions})"


def get_variable_info(db, key):
    """Get the info of a variable as a namedtuple"""
    variable = db.variables[key]
    return DBVariable(variable.datatype, variable.dimensions)


def get_dimension_info(db, key):
    """Get the info of a dimension"""
    dim = db.dimensions[key]
    if dim.isunlimited():
        return 'unlimited'
    else:
        return dim.size

class DatabaseTools(object):
    """Namespace to store general tools to work with the database"""

    @staticmethod
    def get_variables(filename, keys):
        db = load_database(filename)
        result = {key: get_dimension_info(db, key) for key in keys}
        db.close()
        return result

class DatabaseRepresentation(object):
    """Store abstract representation of the Variables and
    Dimensions stored in the Database. This class makes
    also comparisons between different representation 
    straight forward!
    """

    __slots__ = ('variables', 'dimensions', 'unlimited', '_created', '_db', '_handle') 

    def __init__(self, settings):
        self._parse(settings)
        self._created = False
        self._db = None
        self._handle = None

    @classmethod
    def from_db(cls, db):
        """Create DatabaseRepresentation from a database set"""
        variables = {key: get_variable_info(db, key) for key in db.variables.keys()}
        dimensions = {key: get_dimension_info(db, key) for key in db.dimensions.keys()}
        return cls({'variables': variables, 'dimensions': dimensions})

    def __str__(self):
        print("Dimensions: ")
        for item in self['dimensions']:
            print(item)
        for item in self['variables']:
            print(item)
        return ""
    
    def _parse(self, settings):
        """Parse given database"""
        self.dimensions = {}
        self.unlimited = None
        # at least one dimension need to be defined!!!
        for dim_name, dim in settings['dimensions'].items():
            if dim == 'unlimited':
                if self.unlimited is None:
                    self.unlimited = dim
                else:
                    raise Exception("Only a single unlimited dimension allowed!")
            self.dimensions[dim_name] = dim
        # can be that 0 variables are defined, doesnt make sense, but possible
        self.variables = {var_name: variable
                          for var_name, variable in settings.get('variables', {}).items()}

    def __getitem__(self, key):
        if key == 'dimensions':
            return self.dimensions
        elif key == 'variables':
            return self.variables
        else:
            raise KeyError("Only ['dimension', 'variables'] are allowed keys!")
        
    def create_database(self, filename):
        """Create the database from the representation"""
        if self._created is True:
            return self._db, self._handle

        if exists_and_isfile(filename):
            self._db, self._handle = self._load_database(filename)
        else:
            self._db, self._handle = self._init_database(filename)
        # 
        self._created = True
        return self._db, self._handle

    def __eq__(self, rhs):
        """"Compare two representations"""
        assert(isinstance(rhs, self.__class__))
        # check dimensions
        if (set(rhs['dimensions'].keys()) != set(self['dimensions'].keys())):
            return False
        # check that dimensions are exactly the same
        if not all(rhs['dimensions'][dim_name] == self['dimensions'][dim_name] for dim_name in self['dimensions'].keys()):
            return False

        # check variables
        if (set(rhs['variables'].keys()) != set(self['variables'].keys())):
            print("variables keys not the same?")
            return False
        if any(rhs['variables'][var_name] != self['variables'][var_name] for var_name in rhs['variables'].keys()):
            return False
        return True

    def _init_database(self, filename):
        """Create a new database"""
        return create_dataset(filename, self)

    def _load_database(self, filename):
        """Load an existing database and check
           that it is compatable with the existing one"""
        db = load_database(filename)
        ref = DatabaseRepresentation.from_db(db)

        if ref != self:
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
