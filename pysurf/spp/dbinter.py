from abc import abstractmethod
#
import numpy as np

from ..colt import Colt, PluginBase
from ..database.pysurf_db import PySurfDB
from ..database.dbtools import DBVariable
from ..utils.osutils import exists_and_isfile
# logger
from ..logger import get_logger
#
from scipy.linalg import lu_factor, lu_solve


class InterpolatorFactory(PluginBase):
    _is_plugin_factory = True
    _plugins_storage = 'interpolator'


class DataBaseInterpolation(Colt):
    """This class handels all the interaction with the database and
        the interface:
        saves the data and does the interpolation
    """

    _questions = """
        properties = :: list, optional
        write_only = True :: bool
        fit_only = False :: bool
        interpolator = RbfInterpolator :: str
        energy_only = False :: bool
    """
    
    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("interpolator",
                                 {name: interpolator.questions
                                  for name, interpolator in InterpolatorFactory.interpolator.items()})

    def __init__(self, interface, config, natoms, nstates, properties, model=False, logger=None):
        """ """
        if logger is None:
            self.logger = get_logger('db.log', 'database', [])
        else:
            self.logger = logger
        #
        self.write_only = config['write_only']
        self.fit_only = config['fit_only']
        #
        if self.write_only is True and self.fit_only is True:
            raise Exception("Can only write or fit")
        #
        self._interface = interface
        #
        self.natoms = natoms
        self.nstates = nstates
        #
        if config['properties'] is not None:
            properties += config['properties']
        properties += ['crd']
        # setupt database
        self._db = self._create_db(properties, natoms, nstates, model=model)
        self._parameters = get_fitting_size(self._db)
        properties = [prop for prop in properties if prop != 'crd']
        if len(self._db) > 0:
            self.interpolator = InterpolatorFactory.interpolator[config['interpolator'].value](self._db, properties)
        else:
            self.write_only = True


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
            return self.read_last(request)
        if self.write_only is True:
            return self.get_qm(request)
        # do the interpolation
        result, is_trustworthy = self.interpolator.get(request)
        # maybe perform error msg/warning if fitted date is not trustable
        if self.fit_only is True:
            return result
        # do qm calculation
        if is_trustworthy is False:
            return self.get_qm(request)
        return result

    def read_last(self, request):
        for prop in request:
            request.set(prop, self._db.get(prop, -1))
        return request

    def _create_db(self, data, natoms, nstates, filename='db.dat', model=False):
        print('create db model', model)
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


class Interpolator(InterpolatorFactory):

    _register_plugin = False

    def __init__(self, db, properties, savefile=''):
        """important for ShepardInterpolator to set db first!"""
        #
        self.db = db
        #
        if exists_and_isfile(savefile):
            self.interpolators, self.size = self.get_interpolators_from_file(savefile, properties)
        else:
            self.interpolators, self.size = self.get_interpolators(db, properties)

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

    def within_trust_radius(self, crd, radius=0.2):
        is_trustworthy = False
        a = []
        for _crd in self.db['crd']:
            diff = np.linalg.norm((crd - _crd))
            a.append(diff)
            if diff < radius:
                is_trustworthy = True
        return np.asarray(a), is_trustworthy

    def finite_difference_gradient(self, crd, dq=0.01):
        """compute the gradient of the energy  with respect to a crd
           displacement using finite difference method

           TODO: make it usable!
        """
        grad = np.zeros((self.nstates, crd.size), dtype=float)
        #
        crd = crd.resize(3*self.natoms)
        #
        energy = self.interpolators['energy']
        # do loop
        for i in range(crd.size):
            # first add dq
            crd[i] += dq
            en1, _ = energy(crd)
            # first subtract 2*dq
            crd[i] -= 2*dq
            en2, _ = energy(crd)
            # add dq to set crd to origional
            crd[i] += dq
            # compute gradient
            grad[:,i] = (en1 - en2)/(2.0*dq)
        # return gradient
        return grad.resize((self.nstates, self.natoms, 3))


class RbfInterpolator(Interpolator):
    """Basic Rbf interpolator"""

    def get_interpolators(self, db, properties):
        """ """
        A = self._compute_a(db['crd'])
        lu_piv = lu_factor(A)
        return {prop_name: Rbf.from_lu_factors(lu_piv, db[prop_name])
                for prop_name in properties}, len(db)

    def get_interpolators_from_file(self, filename, properties):
        db = Database.load_db(filename)
        out = {}
        for prop_name in db.keys():
            if prop_name == 'size':
                size = db['size']
            if prop_name.endswith('_shape'):
                continue
            if prop_name == 'rbf_epsilon':
                self.epsilon = np.copy(db['rbf_epsilon'])[0]
                continue
            out[prop_name] = Rbf(np.copy(db[prop_name]), tuple(np.copy(db[prop_name+'_shape'])))
        if not all(prop in out for prop in properties):
            raise Exception("Cannot fit all properties")
        return out, size

    def get(self, request):
        """fill request

           Return request and if data is trustworthy or not
        """
        crd, is_trustworthy = self.within_trust_radius(request.crd)
        crd = crd[:self.size]
        crd = weight(crd, self.epsilon)
        for prop in request:
            request.set(prop, self.interpolators[prop](crd))
        #
        return request, is_trustworthy

    def save(self, filename):
        settings = {'dimensions': {}, 'variables': {}}
        dimensions = settings['dimensions']
        variables = settings['variables']

        for prop, rbf in self.interpolators.items():
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
        for prop, rbf in self.interpolators.items():
            db[prop] = rbf.nodes
            db[prop+'_shape'] = rbf.shape
        #
        db['rbf_epsilon'] = self.epsilon

    def _compute_a(self, x):
        size = len(x)
        A = np.empty((size, size), dtype=np.double)
        ndist = 0.0
        # fill the off diagonal
        for i in range(size-1):
            for j in range(i, size):
                res = dist(x[i], x[j])
                A[i, j] = res
                A[j, i] = res
                ndist += res
        # set epsilon
        self.epsilon = ndist/((size-1)*size)
        # fill the diagonal
        for i in range(size):
            A[i, i] = 0.0
        # solve that smarter
        return weight(A, self.epsilon)


class ShepardInterpolator(Interpolator):

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
        return {prop_name: db[prop_name].shape[1:] for prop_name in properties}

    def save(self, filename):
        """Do nothing"""
        print("Warning: ShepardInterpolator does not save files")

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

    def __init__(self, nodes, shape):
        self.nodes = nodes
        self.shape = shape

    @classmethod
    def from_lu_factors(cls, lu_piv, prop):
        nodes, shape = cls._setup(lu_piv, prop)
        return cls(nodes, shape)

    def __call__(self, crd):
        if self.shape == (1,):
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
    return np.sqrt((1.0/epsilon*r)**2 + 1)


def dist(x, y):
    return np.linalg.norm(x-y)
