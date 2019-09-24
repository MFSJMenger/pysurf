"""Tools to store information on the Variables and Dimensions in the Database"""


class DBVariable(object):
    """Store info for Database and easy comparison"""

    __slots__ = ('type', 'dimensions')

    def __init__(self, typ, dim):
        self.type = typ
        self.dimensions = dim

    def __eq__(self, rhs):
        assert(isinstance(rhs, self.__class__))
        if self.type == rhs.type and self.dimensions == rhs.dimensions:
            return True
        return False


def get_variable_info(db, key):
    """Get the info of a variable as a namedtuple"""
    variable = db[key]
    return DBVariable(variable.datatype, variable.dimensions)


def get_dimension_info(db, key):
    """Get the info of a dimension"""
    dim = db[key]
    if dim.isunlimited():
        return 'unlimited'
    else:
        return dim.size


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
        variables = {key: get_variable_info(db.variables, key) for key in db.variables.keys()}
        dimensions = {key: get_dimension_info(db.dimensions, key) for key in db.dimensions.keys()}
        return cls({'variables': variables, 'dimensions': dimensions})
    
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
        self._db, self._handle = create_dataset(filename, settings)
        self._created = True

    def __eq__(self, rhs):
        """"Compare two representations"""
        assert(isinstance(rhs, self.__class__))
        # check dimensions
        if (set(rhs['dimensions'].keys()) != set(self['dimensions'].keys())):
            print('"dimensions" are different between both DatabaseRepresentation!')
            return False
        # check that dimensions are exactly the same
        if not all(rhs['dimensions'][dim_name] == self['dimensions'][dim_name] for dim_name in self['dimensions'].keys()):
            return False

        # check variables
        if (set(rhs['variables'].keys()) != set(self['variables'].keys())):
            return False
        if not all(rhs['variables'][var_name] == self['variables'][var_name] for var_name in rhs['variables'].keys()):
            return False
        return True
