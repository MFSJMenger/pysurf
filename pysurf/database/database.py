from .dbtools import DatabaseRepresentation


class Database(object):
    """"Core Database, can store data given by the DatabaseRepresentation
    
    database is automatically build using a settings dictionary:

    dct = { 'dimensions': {
                 'frame': 'unlimited',
                 'natoms': 3,
                 'three': 3,
                 }, 
             'variables': {
                'coords': DBVariable(np.double, ('frame', 'natoms', 'three')),
             }
          }
    
    """

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
        self._db, self._handle = self._rep.create_database(filename)
        # 
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
