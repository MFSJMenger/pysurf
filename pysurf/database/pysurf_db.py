from collections import namedtuple
import numpy as np

from pysurf.utils import exists_and_isfile
from .dbtools import DatabaseRepresentation, DatabaseGenerator
from .dbtools import load_database as l_db
from .database import Database
from pysurf.system import Molecule, Mode, ModelInfo



class PySurfDB(Database):

    _modes = False
    _molecule = False
    _natoms = False
    _nmodes = False
    _nstates = False
    _atomids = False
    _info = False
    _masses = False
    _crd_equi = False
    _model = False
    _nmodes = False
    _len = False
    _model_info = False

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
        for var in set(cls._variables_molecule.keys()).union(set(cls._variables_model.keys())):
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
            return cls.load_db(filename)
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
        if self._masses is False:
            if 'masses' in self:
                self._masses = np.array(self['masses'])
            else:
                self._masses = None
        return self._masses

    @property
    def modes(self):
        if self._modes is False:
            if 'freqs_equi' in self and 'modes_equi' in self:
                self._modes = [Mode(freq, mode) for freq, mode in zip(np.copy(self['freqs_equi']), np.copy(self['modes_equi']))]
            else:
                self._modes = None
        return self._modes

    @property
    def molecule(self):
        if self.model:
            return None
        if self._molecule is False:
            if 'atomids' in self and 'crd_equi' in self and 'masses' in self:
                self._molecule = Molecule(np.copy(self['atomids']),
                                      np.copy(self['crd_equi']),
                                      np.copy(self['masses']))
            else:
                self._molecule = None
        return self._molecule 

    @property
    def model_info(self):
        if self.model:
            if self._model_info is False:
                if 'crd_equi' in self and 'masses' in self:
                    self._model_info = ModelInfo(np.copy(self['crd_equi']),
                                             np.copy(self['masses']))
                else:
                    self._model_info = None
            return self._model_info
        else:
            return None

    @property
    def system(self):
        if self.model:
            return self.model_info
        else:
            return self.molecule


    @property
    def dimensions(self):
        return self._rep.dimensions

    @property
    def natoms(self):
        if self.model: return None
        if self._natoms is False:
            if 'natoms' in self:
                self._natoms = self.dimensions['natoms']
            else:
                self._natoms = None
        return self._natoms

    @property
    def nstates(self):
        if self._nstates is False:
            if 'nstates' in self.dimensions:
                self._nstates = self.dimensions['nstates']
            else:
                self._nstates = None
        return self._nstates

    @property
    def nmodes(self):
        if self._nmodes is False:
            if 'nmodes' in self.dimensions:
                self._nmodes = self.dimensions['nmodes']
            else:
                self._nmodes = None
        return self._nmodes

    @property
    def len(self):
        if 'crd' in self:
            self._len = len(self['crd'])
        else:
            self._len = None
        return self._len

    def __len__(self):
        return self.len

    @property
    def model(self):
        if self._model is False:
            if 'model' in self:
                self._model =  bool(self['model'][0])
            else:
                self._model = None
        return self._model
