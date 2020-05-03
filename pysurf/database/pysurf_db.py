from collections import namedtuple
import numpy as np

from pysurf.utils import exists_and_isfile
from .dbtools import DatabaseRepresentation, DatabaseGenerator
from .dbtools import load_database as l_db
from .database import Database
from pysurf.molecule import Molecule, Mode



class PySurfDB(Database):

    _modes = None
    _molecule = None
    _natoms = None
    _nstates = None
    _atomids = None
    _info = None
    _masses = None
    _crd_equi = None
    _model = None
    _nmodes = None
    _len = None

    _dimensions_molecule = {
            'frame': 'unlimited',
            'natoms': None,
            'nstates': None,
            'nmodes': None,
            'three': 3,
            'one': 1,
    }
    
    _dimensions_model = {
            'frame': 'unlimited',
            'nstates': None,
            'nmodes': None,
            'one': 1,
    }
    
    _variables_molecule = DatabaseGenerator("""
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
        time      = double :: (frame, one)
        ekin      = double :: (frame, one)
        epot      = double :: (frame, one)
        etot      = double :: (frame, one)
        nacs      = double :: (frame, nstates, nstates, natoms, three)
    """)['variables']

    _variables_model = DatabaseGenerator("""
        [variables]
        crd_equi  = double :: (nmodes)
        freqs_equi= double :: (nmodes)
        modes_equi= double :: (nmodes, nmodes)
        masses    = double :: (nmodes)
        model     = int    :: (one)

        crd       = double :: (frame, nmodes)
        veloc     = double :: (frame, nmodes)
        accel     = double :: (frame, nmodes)
        energy    = double :: (frame, nstates)
        gradient  = double :: (frame, nstates, nmodes)
        fosc      = double :: (frame, nstates)
        currstate = double :: (frame, one)
        time      = double :: (frame, one)
        ekin      = double :: (frame, one)
        epot      = double :: (frame, one)
        etot      = double :: (frame, one)
        nacs      = double :: (frame, nstates, nstates, nmodes)
    """)['variables']


    @classmethod
    def _get_settings(cls, variables, model):
        out = {'variables': {}, 'dimensions': {}}
        varis = out['variables']
        dims = out['dimensions']
        
        if model: 
            refvar = cls._variables_model
            refdim = cls._dimensions_model
        else: 
            refvar = cls._variables_molecule
            refdim = cls._dimensions_molecule

        for var in variables:
            try:
                varis[var] = refvar[var]
            except:
                raise Exception(f"Variable {var} unknown")
            for dim in varis[var].dimensions:
                try:
                    dims[dim] = refdim[dim]
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
        settings = cls._get_settings(data, model)
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

    def add_reference_entry(self, system, modes, model):
        self.set('model', int(model))
        if model is False:
            self.set('atomids', system.atomids)
        self.set('masses', system.masses)
        self.set('modes_equi', np.array([mode.displacements for mode in modes]))
        self.set('freqs_equi', np.array([mode.freq for mode in modes]))
        # add equilibrium values
        self.set('crd_equi', system.crd)

    @property
    def masses(self):
        if self._masses is None:
            self._masses = np.array(self['masses'])
        return self._masses

    @property
    def modes(self):
        if self._modes is None:
            self._modes = [Mode(freq, mode) for freq, mode in zip(np.copy(self['freqs_equi']), np.copy(self['modes_equi']))]
        return self._modes

    @property
    def molecule(self):
        if self.model:
            return None
        if self._molecule is None:
            self._molecule = Molecule(np.copy(self['atomids']),
                                      np.copy(self['crd_equi']),
                                      np.copy(self['masses']))
        return self._molecule 

    @property
    def dimensions(self):
        return self._rep.dimensions

    @property
    def natoms(self):
        if self._natoms is None:
            self._natoms = self.dimensions['natoms']
        return self._natoms

    @property
    def nstates(self):
        if self._nstates is None:
            self._nstates = self.dimensions['nstates']
        return self._nstates

    @property
    def len(self):
        if self._len is None:
            self._len = len(self['crd'])
        return self._len

    @property
    def model(self):
        if self._model is None:
            self._model =  bool(self['model'][0])
        return self._model
