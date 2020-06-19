from pysurf import Interpolator


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
