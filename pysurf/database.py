from netCDF4 import Dataset

from .utils.osutils import exists_and_isfile


class Database(object):

    def __init__(self, filename, settings=None):
        """Initialize new Database,
        if db exisits:
           load existing database
           check that the settings of the old database 
           are the same with the once used in the loading
           routine
        else:
           create new database
        """
        if exists_and_isfile(filename):
            self._load_database(filename, settings)
        else:
            self._init_database(filename, settings)

    def _init_database(self, filename, settings):
        """Create a new database"""
        pass

    def _load_database(self, filename, settings):
        """Load an existing database and check 
           that it is compatable with the existing one"""
        pass
