from ..colt import Colt
from ..database.database import Database
# logger
from ..logger import get_logger


class DataBaseInterpolation(Colt):
    """This class handels all the interaction with the database and
        the interface:
        saves the data and does the interpolation
    """

    _questions = """
        properties = none :: list
        write_only = True :: bool
    """

    def __init__(self, interface, config, natoms, nstates, dynamics_properties, logger=None):
        """

        """
        if logger is None:
            self.logger = get_logger('db.log', 'database', [])
        else:
            self.logger = logger
        #
        self._write_only = config['write_only']
        #
        self._interface = interface
        self.properties = dynamics_properties
        self.natoms = natoms
        self.nstates = nstates
        #
        self._db = self._create_db()

    def _create_db(self, filename='db.dat'):
        data = f"""
        [dimensions]
        frame = unlimited
        natoms = {self.natoms}
        nstates = {self.nstates}
        three = 3

        [variables]
        crd = double :: (frame, natoms, three)
        energy = double :: (frame, nstates)
        gradient = double :: (frame, nstates, natoms, three)
        """
        return Database(filename, data)

    def get(self, request):
        """Get result of request"""
        #
        result = self._interface.get(request)
        #
        self._db.append('crd', result['crd'])
        self._db.append('energy', result['energy'])
        self._db.append('gradient', result['gradient'].data)
        #self._db.append('dipol', result['dipol'])
        self._db.increase
        result = {key: result[key] for key in request.keys() if key != 'crd'}
        return result
