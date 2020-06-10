from abc import abstractmethod
from functools import reduce
#
import numpy as np

from ..colt import Colt, PluginBase
from ..database.pysurf_db import PySurfDB
from ..database.database import Database
from ..database.dbtools import DBVariable
from ..utils.osutils import exists_and_isfile
# logger
from ..logger import get_logger
#
from scipy.linalg import lu_factor, lu_solve
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.spatial import cKDTree

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures


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



class RegInterpolator(Interpolator):
    """Basic Rbf interpolator"""

    _questions = """
        trust_radius_general = 0.75 :: float
        trust_radius_ci = 0.25 :: float
        energy_threshold = 0.02 :: float
        regression_order = 2 :: int
    """

    @classmethod
    def from_config(cls, config, db, properties, logger, energy_only, weightsfile, crdmode, fit_only):
        trust_radius_general = config['trust_radius_general']
        trust_radius_CI = config['trust_radius_ci']
        energy_threshold = config['energy_threshold']
        order = config['regression_order']
        #
        return cls(db, properties, logger, energy_only=energy_only, weightsfile=weightsfile,
                   crdmode=crdmode, trust_radius_general=trust_radius_general,
                   trust_radius_CI=trust_radius_CI, energy_threshold=energy_threshold, order=order, fit_only=fit_only)

    def __init__(self, db, properties, logger, energy_only=False, weightsfile=None, crdmode='cartesian',
                 trust_radius_general=0.75, trust_radius_CI=0.25, energy_threshold=0.02, order=2, fit_only=False):

        self.trust_radius_general = trust_radius_general
        self.trust_radius_CI = trust_radius_CI
        self.energy_threshold = energy_threshold
        self.order = order
        super().__init__(db, properties, logger, energy_only, weightsfile, crdmode=crdmode, fit_only=fit_only)


    def get_interpolators(self, db, properties):
        return {prop_name: Regression(self.crds, db[prop_name], self.order)
                for prop_name in properties}, len(db['crd'])

    def get_interpolators_from_file(self, filename, properties):
        pass

    def get(self, request):
        """fill request

           Return request and if data is trustworthy or not
        """
        if self.crdmode == 'internal':
            crd = internal(request.crd)
        else:
            crd = request.crd
        #
        _, trustworthy = within_trust_radius(crd, self.crds, radius=self.trust_radius_general, radius_ci=self.trust_radius_CI)
#       crd = crd[:self.size]


        for prop in request:
            request.set(prop, self.interpolators[prop](crd, request))
        #
        diffmin = np.min(np.diff(request['energy']))
        #compare energy differences with threshold from user
        if diffmin < self.energy_threshold:
            self.logger.info(f"Small energy gap of {diffmin}. Within CI radius: " + str(trustworthy[1]))
            is_trustworthy = trustworthy[1]
        else:
            self.logger.info('Large energy diffs. Within general radius: ' + str(trustworthy[0]))
            is_trustworthy = trustworthy[0]
        return request, is_trustworthy


    def save(self, filename):
        pass

    def _train(self):
        self.crds = self.get_crd()
        #
        for name, interpolator in self.interpolators.items():
            if isinstance(interpolator, Regression):
                interpolator.train()
    
    def loadweights(self, filename):
        pass


class Regression:
    def __init__(self, crds, values, order):
#       self.crds = np.copy(crds)
        self.crds = crds
        self.shape_crds = self.crds.shape
        self.crds.resize(self.shape_crds[0], int(self.crds.size/self.shape_crds[0]))
        self.values = np.copy(values)
        self.shape_values = values.shape
        self.values.resize(self.shape_values[0], int(self.values.size/self.shape_values[0]))
        self.poly = PolynomialFeatures(degree=order)

    def train(self):
        self.crds_poly = self.poly.fit_transform(self.crds)
        self.regr = LinearRegression()
        self.regr.fit(self.crds_poly, self.values)

    def __call__(self, crd, request):
        crd = np.copy(crd)
        crd = self.poly.fit_transform(crd.reshape((1,crd.size)))
        res = self.regr.predict(crd)
        res.resize(self.shape_values[1:])
        return res


class RbfInterpolator(Interpolator):
    """Basic Rbf interpolator"""

    _questions = """
        trust_radius_general = 0.75 :: float
        trust_radius_ci = 0.25 :: float
        energy_threshold = 0.02 :: float
        epsilon = :: float, optional
    """

    @classmethod
    def from_config(cls, config, db, properties, logger, energy_only, weightsfile, crdmode, fit_only):
        trust_radius_general = config['trust_radius_general']
        trust_radius_CI = config['trust_radius_ci']
        energy_threshold = config['energy_threshold']
        epsilon = config['epsilon']
        #
        return cls(db, properties, logger, energy_only=energy_only, weightsfile=weightsfile,
                   crdmode=crdmode, trust_radius_general=trust_radius_general,
                   trust_radius_CI=trust_radius_CI, energy_threshold=energy_threshold, fit_only=fit_only, epsilon=epsilon)

    def __init__(self, db, properties, logger, energy_only=False, weightsfile=None, crdmode='cartesian', fit_only=False,
            trust_radius_general=0.75, trust_radius_CI=0.25, energy_threshold=0.02, epsilon=None):

        self.trust_radius_general = trust_radius_general
        self.trust_radius_CI = trust_radius_CI
        self.energy_threshold = energy_threshold
        self.trust_radius = (self.trust_radius_general + self.trust_radius_CI)/2.
        if epsilon is not None:
            self.epsilon = epsilon
        else:
            self.epsilon = trust_radius_CI
        super().__init__(db, properties, logger, energy_only, weightsfile, crdmode=crdmode, fit_only=fit_only)


    def get_interpolators(self, db, properties):
        """ """
        A = self._compute_a(self.crds)
        lu_piv = lu_factor(A)
        return {prop_name: Rbf.from_lu_factors(lu_piv, db[prop_name], self)
                for prop_name in properties}, len(db)

    def get_interpolators_from_file(self, filename, properties):
        db = Database.load_db(filename)
        out = {}
        for prop_name in db.keys():
            if prop_name.endswith('_shape'):
                continue
            if prop_name == 'rbf_epsilon':
                self.epsilon = np.copy(db['rbf_epsilon'])[0]
                continue
            out[prop_name] = Rbf(np.copy(db[prop_name]), tuple(np.copy(db[prop_name+'_shape'])), self)
        if not all(prop in out for prop in properties):
            raise Exception("Cannot fit all properties")
        return out

    def get(self, request):
        """fill request

           Return request and if data is trustworthy or not
        """
        if self.crdmode == 'internal':
            crd = internal(request.crd)
        else:
            crd = request.crd
        #
        _, trustworthy = within_trust_radius(crd, self.crds, radius=self.trust_radius_general, radius_ci=self.trust_radius_CI)

        for prop in request:
            request.set(prop, self.interpolators[prop](crd, request))
        #
        diffmin = np.min(np.diff(request['energy']))
        #compare energy differences with threshold from user
        if diffmin < self.energy_threshold:
            self.logger.info(f"Small energy gap of {diffmin}. Within CI radius: " + str(trustworthy[1]))
            is_trustworthy = trustworthy[1]
        else:
            self.logger.info('Large energy diffs. Within general radius: ' + str(trustworthy[0]))
            is_trustworthy = trustworthy[0]
        return request, is_trustworthy

    def loadweights(self, filename):
        """Load existing weights"""
        db = Database.load_db(filename)
        for prop, rbf in self.interpolators.items():
            if prop == 'gradient' and self.energy_only is True:
                continue
            if prop not in db:
                raise Exception("property needs to be implemented")
            rbf.nodes = np.copy(db[prop])
            rbf.shape = np.copy(db[prop+'_shape'])
        self.epsilon = np.copy(db['rbf_epsilon'])

    def save(self, filename):
        settings = {'dimensions': {}, 'variables': {}}
        dimensions = settings['dimensions']
        variables = settings['variables']

        for prop, rbf in self.interpolators.items():
            if not isinstance(rbf, Rbf):
                continue
            lshape = len(rbf.shape)
            dimensions[str(lshape)] = lshape
            for num in rbf.nodes.shape:
                dimensions[str(num)] = num
            variables[prop] = DBVariable(np.double, tuple(str(num) for num in rbf.nodes.shape))
            variables[prop+'_shape'] = DBVariable(np.int, tuple(str(lshape)))
        dimensions['1'] = 1
        variables['rbf_epsilon'] = DBVariable(np.double, ('1',))
        #
        db = Database(filename, settings)
        #
        for prop, rbf in self.interpolators.items():
            if not isinstance(rbf, Rbf):
                continue
            db[prop] = rbf.nodes
            db[prop+'_shape'] = rbf.shape
        #
        db['rbf_epsilon'] = self.epsilon

    def _train(self):
        """set rbf weights, based on the current crds"""
#       self.crds = self.get_crd()
        A = self._compute_a(self.crds)
        lu_piv = lu_factor(A)
        #
        for name, interpolator in self.interpolators.items():
            if isinstance(interpolator, Rbf):
                interpolator.update(lu_piv, self.db[name])

    def _compute_a(self, x):
        #
        shape = x.shape
        if len(shape) == 3:
            dist = pdist(x.reshape((shape[0], shape[1]*shape[2])))
        else:
            dist = pdist(x)
        A = squareform(dist)
        return weight(A, self.epsilon)


def dim_norm(crd1, crd2):
    return np.max(np.abs(crd1-crd2))


class ShepardInterpolator(Interpolator):

    def __init__(self, db, properties, logger, energy_only=False, weightsfile=None, crdmode='cartesian', fit_only=False):
        super().__init__(db, properties, logger, energy_only, weightsfile, crdmode=crdmode, fit_only=fit_only)
        self.crds = self.get_crd()

    def get(self, request):
        """fill request and return True
        """
        #
        weights, is_trustworthy = self._get_weights(request.crd)
        # no entries in db...
        if weights is None:
            return request, False
        #
        for prop in request:
            request.set(prop, self._get_property(weights, prop))
        #
        return request, is_trustworthy

    def get_interpolators(self, db, properties):
        return {prop_name: db[prop_name].shape[1:] for prop_name in properties}, len(db)

    def save(self, filename):
        """Do nothing"""

    def loadweights(self, filename):
        """Do nothing"""

    def _train(self):
        pass

    def get_interpolators_from_file(self, filename, properties):
        return {prop_name: self.db[prop_name].shape[1:] for prop_name in properties}

    def _get_property(self, weights, prop):
        entries = db[prop]
        shape = self.interpolators[prop]
        res = np.zeros(shape, dtype=np.double)
        for i, value in enumerate(entries):
            res += weights[i]*value
        res = res/np.sum(weights)
        if shape == (1,):
            return res[0]
        return res

    def _get_weights(self, crd, trust_radius=0.2):
        """How to handle zero division error"""
        exact_agreement = False
        crds = db['crd']
        size = len(crds)
        if size == 0:
            return None, False
        #
        weights = np.zeros(size, dtype=np.double)
        #
        is_trustworthy = False
        for i in range(size):
            diff = np.linalg.norm((crd-crds[i]))**2
            if diff < trust_radius:
                is_trustworthy = True
            if round(diff, 6) == 0:
                exact_agreement = i
            else:
                weights[i] = 1./diff
        if exact_agreement is False:
            return weights, is_trustworthy
        #
        weights.fill(0.0)
        weights[exact_agreement] = 1.0
        return weights, is_trustworthy


class Rbf:

    def __init__(self, nodes, shape, parent):
        self.nodes = nodes
        self.shape = shape
        self.parent = parent

    def update(self, lu_piv, prop):
        self.nodes, self.shape = self._setup(lu_piv, prop)

    @classmethod
    def from_lu_factors(cls, lu_piv, prop, parent):
        nodes, shape = cls._setup(lu_piv, prop)
        return cls(nodes, shape, parent)

    def __call__(self, crd, request):
        shape = self.parent.crds.shape
        if len(shape) == 3:
            dist = cdist([np.array(crd).flatten()], self.parent.crds.reshape((shape[0], shape[1]*shape[2])))
        else:
            dist = cdist([np.array(crd).flatten()], self.parent.crds)
        crd = weight(dist, self.parent.epsilon)
        if len(self.shape) == 1 and self.shape[0] == 1:
            return np.dot(crd, self.nodes).reshape(self.shape)[0]
        return np.dot(crd, self.nodes).reshape(self.shape)

    @staticmethod
    def _setup(lu_piv, prop):
        prop = np.array(prop)
        shape = prop.shape
        size = shape[0]
        dim = 1
        for i in shape[1:]:
            dim *= i
        #
        prop = prop.reshape((size, dim))
        #
        if dim == 1:
            nodes = lu_solve(lu_piv, prop)
        else:
            nodes = np.zeros((size, dim), dtype=prop.dtype)
            for i in range(dim):
                nodes[:,i] = lu_solve(lu_piv, prop[:,i])
        return nodes, shape[1:]


def weight(r, epsilon):
#    return r
    return np.sqrt((1.0/epsilon*r)**2 + 1)


def dist(x, y):
    return np.linalg.norm(x-y)
