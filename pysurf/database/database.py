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

    __slots__ = ('filename', '_rep', '_db', '_handle', 'closed', '_icurrent')

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
        # 
        self._icurrent = None

    def __getitem__(self, key):
        return self._handle.get(key, None)

    def __contains__(self, key):
        return key in self._handle

    def get(self, key, ivalue):
        variable = self._handle[key]
        if variable.shape[0] > ivalue:
            return variable[ivalue]

    def set(self, ivalue, key, value):
        if key in self:
            self._handle[key][ivalue, :] = value

    def get_dimension_size(self, key):
        if key in self._db.dimensions.keys():
            return self._db.dimensions[key].size

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

    def close(self):
        if self.closed is False:
            self._db.close()
            self.closed = True

    def __del__(self):
        self.close()
