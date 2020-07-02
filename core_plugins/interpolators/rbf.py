import numpy as np
#
from scipy.linalg import lu_factor, lu_solve
from scipy.spatial.distance import cdist, pdist, squareform
#
from pysurf import Interpolator
from pysurf.spp import within_trust_radius, internal
from pysurf.database.dbtools import DBVariable
from pysurf.database.database import Database
#
from codetiming import Timer

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

#    @Timer(name="get")
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

#    @Timer(name="train")
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
