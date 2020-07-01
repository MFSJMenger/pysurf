from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from scipy.spatial.distance import cdist
#
from pysurf import Interpolator
from pysurf.spp import within_trust_radius, internal


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
