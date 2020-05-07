from ..colt import Colt
from ..database.pysurf_db import PySurfDB
# logger
from ..logger import get_logger


class DataBaseInterpolation(Colt):
    """This class handels all the interaction with the database and
        the interface:
        saves the data and does the interpolation
    """

    _questions = """
        properties = :: list, optional
        write_only = True :: bool
        fit_only = False :: bool
        interpolation = shepard :: str :: [shepard]
        energy_only = False :: bool
    """

    def __init__(self, interface, config, natoms, nstates, properties, logger=None):
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
        # setupt database
        self._db = self._create_db(properties, natoms, nstates)
        self._parameters = get_fitting_size(self._db)
        #
        self.interpolator = Interpolator(self._db, self._parameters, 
                                         natoms, nstates, 
                                         config, self.logger)

    def get_qm(self, request):
        """Get result of request and append it to the database"""
        #
        result = self._interface.get(request)
        #
        for prop in self._parameters:
            self._db.append(prop, result['prop'])
        #
        self._db.increase
        return result

    def get(self, request):
        """answer request"""
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

    def _create_db(self, data, natoms, nstates, filename='db.dat', is_model=False):
        if is_model is False:
            return PySurfDB.generate_database(filename, data=data, dimensions={'natoms': natoms, 'nstates': nstates})
        return PySurfDB.generate_database(filename, data=data, dimensions={'nmodes': natoms, 'nstates': nstates})


class Interpolator(Colt):
    """Factory class for all Interpolation"""

    interpolator_schemes = {'shepard': ShepardInterpolator}

    def __init__(self, db, parameters, natoms, nstates, config, logger):
        self.db = db
        self.logger = logger
        self.nstates = nstates
        self.natoms = natoms
        # The interpolate gradient key has been added in DBInter
        self.energy_only = config['energy_only']

        self.logger.info(f"set up interpolation for {config['interpolation']}-Interpolationscheme")
        interpolation_scheme = self.interpolator_schemes[config['interpolation']]
        self.interpolators = self._get_interpolators(parameters, interpolation_scheme, db)
        self.logger.info('successfully set up interpolation')

    def get(self, request):
        """Answers request using the interpolators
        
           Returns: request, bool
           bool = True: data is within trust radius
                  False: data is outside trust radius
        """
        crd = request['crd']
        for prop in request.keys():
            if prop == 'crd':
                continue
            if prop == 'gradient' and self.energy_only is True:
                request[prop] = self.finite_difference_gradient(crd)
                continue
            request[prop] = self.interpolators[prop](crd)
        return request, True
        
    def finite_difference_gradient(self, crd, dq=0.01):
        """compute the gradient of the energy  with respect to a crd 
           displacement using finite difference method
        """
        grad = np.zeros((self.nstates, crd.size), dtype=float)
        #
        crd = crd.resize(3*self.natoms)
        energy = self.interpolators['energy']
        # do loop
        for i in range(crd.size):
            # first add dq
            crd[i] += dq
            en1 = energy(crd)
            # first subtract 2*dq
            crd[i] -= 2*dq
            en2 = energy(crd)
            # add dq to set crd to origional
            crd[i] += dq
            # compute gradient
            grad[:,i] = (en1 - en2)/(2.0*dq)
        # return gradient
        return grad.resize((self.nstates, self.natoms, 3))

    def _get_interpolators(self, parameters, scheme, db):
        return {param: scheme(db['crd'], db[param], size) for param, size in parameters.items() if param != 'crd'}


class ShepardInterpolator:

    def __init__(self, crds, values, size):
        self.crds = crds
        self.values = values
        self.size = size

    def __call__(self, crd):
        weights = self._get_weights(crd)
        res = np.zeros(self.size, dtype=np.double)
        for i, value in enumerate(self.values):
            res += weights[i]*value
        res = res/np.sum(weights)
        return res

    def _get_weights(self, crd):
        """How to handle zero division error"""
        size = len(self.crds)
        #
        weights = np.zeros(size, dtype=np.double)
        #
        for i in range(size):
            diff = np.linalg.norm((crd-self.crds[i]))**2
            if round(diff, 6) == 0:
                exact_agreement = i
            else:
                weights[i] = 1./diff
        if exact_agreement is False:
            return weights
        #
        weights.fill(0.0)
        weights[exact_agreement] = 1.0
        return weights


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
