from abc import abstractmethod
#
from scipy.spatial.distance import cdist, pdist
import numpy as np

from ..colt import Colt, PluginBase
from ..database.pysurf_db import PySurfDB
from ..utils.osutils import exists_and_isfile
# logger
from ..logger import get_logger


def internal(crd):
    return pdist(crd)
#    return np.array([1.0/ele for ele in pdist(crd)])


def internal_coordinates(crds):
    return np.array([internal(crd) for crd in crds])


class DataBaseInterpolation(Colt):
    """This class handels all the interaction with the database and
        the interface:
        saves the data and does the interpolation
    """

    _questions = """
        # additional properties to be fitted
        properties = :: list, optional
        # only write
        write_only = True :: bool
        # name of the database
        database = db.dat :: file
    """

    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_block("interpolator", Interpolator.questions)

    @classmethod
    def from_config(cls, config, interface, natoms, nstates, properties, model=False, logger=None):
        return cls(interface, config, natoms, nstates, properties, model=model, logger=logger)

    def __init__(self, interface, config, natoms, nstates, properties, model=False, logger=None):
        """ """
        self.config = config
        if logger is None:
            self.logger = get_logger('db.log', 'database', [])
        else:
            self.logger = logger
        #
        self.write_only = config['write_only']

        #
        self._interface = interface
        #
        self.natoms = natoms
        self.nstates = nstates
        #
        if config['properties'] is not None:
            properties += config['properties']
        properties += ['crd']
        # setup database
        self._db = self._create_db(properties, natoms, nstates, model=model, filename=config['database'])
        self._parameters = get_fitting_size(self._db)
        properties = [prop for prop in properties if prop != 'crd']
        self.properties = properties
        if config['write_only'] is False:
            self.interpolator = Interpolator.setup_from_config(config['interpolator'], self._db,
                                                    properties,
                                                    logger=self.logger)
            self.fit_only = self.interpolator.fit_only
            #
            if self.write_only is True and self.fit_only is True:
                raise Exception("Can only write or fit")
        else:
            self.write_only = True
            self.fit_only = False


    def get_qm(self, request):
        """Get result of request and append it to the database"""
        #
        result = self._interface.get(request)
        #
        for prop, value in result.iter_data():
            self._db.append(prop, value)
        self._db.append('crd', result.crd)
        #
        self._db.increase
        return result

    def get(self, request):
        """answer request"""
        if request.same_crd is True:
            return self.old_request
        self.old_request = self._get(request)
        return self.old_request

    def _get(self, request):
        """answer request"""
        if self.write_only is True:
            return self.get_qm(request)
        # do the interpolation
        result, is_trustworthy = self.interpolator.get(request)
        # maybe perform error msg/warning if fitted date is not trustable
        if self.fit_only is True:
            if is_trustworthy is False:
                self.logger.warning('Interpolated result not trustworthy, but used as fit_only is True')
            return result
        # do qm calculation
        if is_trustworthy is False:
            self.logger.info('Interpolated result is not trustworthy and QM calculation is started')
            return self.get_qm(request)
        self.logger.info('Interpolated result is trustworthy and returned')
        return result

    def read_last(self, request):
        for prop in request:
            request.set(prop, self._db.get(prop, -1))
        return request

    def _create_db(self, data, natoms, nstates, filename='db.dat', model=False):
        if model is False:
            return PySurfDB.generate_database(filename, data=data, dimensions={'natoms': natoms, 'nstates': nstates, 'nactive': nstates}, model=model)
        return PySurfDB.generate_database(filename, data=data, dimensions={'nmodes': natoms, 'nstates': nstates, 'nactive': nstates}, model=model)


def get_fitting_size(db):
    """We only fit unlimeted data"""
    out = {}
    for variable in db.get_keys():
        ndim = 1
        dims = db.get_dimension(variable)
        if not dims[0].isunlimited():
            continue
        for dim in dims[1:]:
            ndim *= dim.size
        out[variable] = ndim
    return out


class InterpolatorFactory(PluginBase):
    _is_plugin_factory = True
    _plugins_storage = 'interpolator'


class Interpolator(InterpolatorFactory):

    _questions="""
    weights_file = :: file, optional
    # only fit to the existing data
    fit_only = False :: bool
    # select interpolator
    interpolator = RbfInterpolator :: str
    # if true: compute gradient numerically
    energy_only = False :: bool
    crdmode = internal :: str :: [internal, cartesian]
    """

    _register_plugin = False

    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_cases("interpolator",
                                 {name: interpolator.questions
                                  for name, interpolator in cls.plugins.items()})
    
    @classmethod
    def setup_from_config(cls, config, db, properties, logger):
        return InterpolatorFactory.plugin_from_config(config['interpolator'], db,
                                                    properties,
                                                    logger=logger,
                                                    energy_only=config['energy_only'],
                                                    weightsfile=config['weights_file'],
                                                    crdmode=config['crdmode'],
                                                    fit_only=config['fit_only'])
    
    @classmethod
    def from_config(cls, config, db, properties, logger, energy_only, weightsfile, crdmode, fit_only):
        return cls(db, properties, logger, energy_only, weightsfile, crdmode, fit_only)

    def __init__(self, db, properties, logger, energy_only=False, weightsfile=None, crdmode=False, fit_only=False):
        """important for ShepardInterpolator to set db first!"""
        #
        self.crds = None
        self.logger = logger
        self.db = db
        self.nstates = self.db.get_dimension_size('nstates')
        self.energy_only = energy_only
        self.fit_only = fit_only
        self.weightsfile = weightsfile
        #
        self.crdmode = crdmode
        self.crds = self.get_crd()
        #
        if energy_only is True:
            properties = [prop for prop in properties if prop != 'gradient']
        #
        if exists_and_isfile(weightsfile):
            print('weightsfile', weightsfile)
            self.interpolators = self.get_interpolators_from_file(weightsfile, properties)
        else:
            self.interpolators, self.size = self.get_interpolators(db, properties)
        #
        if energy_only is True:
            self.interpolators['gradient'] = self.finite_difference_gradient
        # train the interpolator!
        self.train(weightsfile)

    def get_crd(self):
        if self.crdmode == 'internal':
            crds = internal_coordinates(np.copy(self.db['crd']))
        else:
            crds = np.copy(self.db['crd'])
        return crds

    @abstractmethod
    def get(self, request):
        """fill request

           Return request and if data is trustworthy or not
        """

    @abstractmethod
    def get_interpolators(self, db, properties):
        """ """

    @abstractmethod
    def save(self, filename):
        """Save weights"""

    @abstractmethod
    def get_interpolators_from_file(self, filename, properties):
        """setup interpolators from file"""

    @abstractmethod
    def _train(self):
        """train the interpolators using the existing data"""

    @abstractmethod
    def loadweights(self, filename):
        """load weights from file"""

    def train(self, filename=None, always=False):
        if filename == '':
            filename = None
        # train normally
        if exists_and_isfile(filename):
            self.loadweights(filename)
            return
        else:
            self._train()
        # save weights
        if filename is not None:
            self.save(filename)

    def update_weights(self):
        """update weights of the interpolator"""
        self.train()

    def finite_difference_gradient(self, crd, request, dq=0.01):
        """compute the gradient of the energy  with respect to a crd
           displacement using finite difference method
        """
        crd = request.crd
        grad = np.zeros((self.nstates, crd.size), dtype=float)
        #
        shape = crd.shape
        crd.resize(crd.size)
        #
        energy = self.interpolators['energy']
        # do loop
        for i in range(crd.size):
            # first add dq
            crd[i] += dq
            if self.crdmode == 'internal':
                crd_here = internal(crd.reshape(shape))
            else:
                crd_here = crd
            en1 = energy(crd_here, request)
            # first subtract 2*dq
            crd[i] -= 2*dq
            if self.crdmode == 'internal':
                crd_here = internal(crd.reshape(shape))
            else:
                crd_here = crd
            en2 = energy(crd_here, request)
            # add dq to set crd to origional
            crd[i] += dq
            # compute gradient
            grad[:,i] = (en1 - en2)/(2.0*dq)
        # return gradient
        crd.resize(shape)
        grad.resize((self.nstates, *shape))
        return grad


def within_trust_radius(crd, crds, radius, metric='euclidean', radius_ci=None):
    is_trustworthy_general = False
    is_trustworthy_CI = False
    shape = crds.shape
    crd_shape = crd.shape
    crd.resize((1, crd.size))
    if len(shape) == 3:
        crds.resize((shape[0], shape[1]*shape[2]))
    dist = cdist(crd, crds, metric=metric)
    crds.resize(shape)
    crd.resize(crd_shape)
    if np.min(dist) < radius:
        is_trustworthy_general = True
    if radius_ci is not None:
        if np.min(dist) < radius_ci:
            is_trustworthy_CI = True
        return dist[0], (is_trustworthy_general, is_trustworthy_CI)
    else:
        return dist[0], is_trustworthy_general
