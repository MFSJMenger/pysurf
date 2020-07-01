import numpy as np
#
from scipy.spatial import cKDTree
#
from pysurf import Interpolator
from pysurf.spp import internal


class NearestNeighborInterpolator(Interpolator):
    """Nearest Neighbor Interpolator"""

    _questions = """
        trust_radius_general = 0.75 :: float
        trust_radius_ci = 0.25 :: float
        energy_threshold = 0.02 :: float
        norm = euclidean :: str :: [manhattan, euclidean, max]
    """

    @classmethod
    def from_config(cls, config, db, properties, logger, energy_only, weightsfile, crdmode, fit_only):
        trust_radius_general = config['trust_radius_general']
        trust_radius_CI = config['trust_radius_ci']
        energy_threshold = config['energy_threshold']
        #
        # convert input for norm in corresponding input (p-Norm) for the cKDTree
        # for more information go to the cKDTree.query documentation
        if config['norm'] == 'manhattan':
            norm = 1
        elif config['norm'] == 'max':
            norm = 'infinity'
        else:
            norm = 2
        #
        return cls(db, properties, logger, energy_only=energy_only, weightsfile=weightsfile,
                   crdmode=crdmode, trust_radius_general=trust_radius_general,
                   trust_radius_CI=trust_radius_CI, energy_threshold=energy_threshold,
                   fit_only=fit_only, norm=norm)

    def __init__(self, db, properties, logger, energy_only=False, weightsfile=None,
                 crdmode='cartesian', fit_only=False, trust_radius_general=0.75,
                 trust_radius_CI=0.25, energy_threshold=0.02, norm='euclidean'):
        """ Storing user input as internal variables and call the init method of the 
            interpolator factory

            Parameters
            ----------

                db:
                    databse containing the datasets, on which the interpolation is
                    based on

                properties: list
                    properties (e.g. ['energy', 'gradient']) that should be fitted 

                logger:
                     logger to log any incident

                energy_only: bool, optional
                    if energy_only is True, gradients are derived from the energy surface

                weightsfile: str, optional
                    filepath, where to save the weights. Not used in the case of the 
                    NearestNeighborInterpolator, but needed for the overall framework.

                crdmode: str, optional
                    Variable to determine whether a coordinate transformation is applied before 
                    fitting.

                fit_only: bool, optional
                    Flag to determine, whether no new QM calculations are performed

                trust_radius_general: float, optional
                    radius to determine whether fitted result is trustworthy in regions of a 
                    large energy gap

                trust_radius_CI: float, optional
                    radius to determine whether fitted result is trustworthy in regions of a 
                    small energy gap

                energy_threshold: float, optional
                    Threshold to distinguish regions of small and large energy gaps.

                norm: str, optional
                    Determining the norm for the nearest neighbor search. 'manhattan' corresponds
                    to the 1-norm, 'euclidean' is the 2-norm, and 'max' is the infinity norm. 
        """
        #
        self.trust_radius_general = trust_radius_general
        self.trust_radius_CI = trust_radius_CI
        self.energy_threshold = energy_threshold
        self.tree = None
        self.norm = norm
        #
        super().__init__(db, properties, logger, energy_only, weightsfile, crdmode=crdmode,
                         fit_only=fit_only)

    def get_interpolators(self, db, properties):
        """ For each property a separate interpolator is set up to be consistent with the
            PySurf Interpolator Framework. To avoid setting up several cKDTrees the same 
            tree is used for all interpolators.

            Parameters:
            -----------
                db:
                    databse containing the datasets, on which the interpolation is
                    based on

                properties: list
                    properties (e.g. ['energy', 'gradient']) that should be fitted 

            
            Returns
            --------
                dictionary containing the property name as key and the corresponding
                interpolator as value.
        """
        #
        if self.tree is None:
            self.tree = cKDTree(self.crds)
        return {prop_name: NNInterpolator(db, self.tree, prop_name, norm=self.norm)
                for prop_name in properties}, len(db)


    def get_interpolators_from_file(self, filename, properties):
        """ Specifically for Machine Learning algorithms interpolators can be loaded from a file.
            However, in the case of NearestNeighborInterpolation that is not used and implemented.
            To be consistent, the get_interpolators method is called
            
            Parameters:
            -----------
                filename: string
                    filepath where information of interpolators is stored. Not used here!

                properties: list
                    properties (e.g. ['energy', 'gradient']) that should be fitted 
            

            Returns
            --------
                dictionary containing the property name as key and the corresponding
                interpolator as value.
        """

        self.logger.warning("NearestNeighborInterpolator cannot be started from a file. " +
                            "Interpolators are set up from database")
        return self.get_interpolators(self.db, properties)

#    @Timer(name="get")
    def get(self, request):
        """ Fill request and return request and if data is trustworthy or not

            Parameters:
            -----------
                request:
                    request instance of request class, which is the standardized communication
                    between the spp and its clients
            
            Returns:
            -----------
                request:
                    The same request instance as the input parameter, but the desired information
                    is filled in.
        """
        #
        # Convert coordinate into desired format
        if self.crdmode == 'internal':
            crd = internal(request.crd)
        else:
            crd = request.crd
        #
        # Make nearest neighbor search once and pass it to all interpolators
        dist, idx = self.tree.query(crd, p=self.norm)
        for prop in request:
            request.set(prop, self.interpolators[prop](crd, request, idx))
        #
        # Determine whether result is trustworthy, using the trust radii
        diffmin = np.min(np.diff(request['energy']))
        is_trustworthy = False
        if diffmin < self.energy_threshold:
            if dist < self.trust_radius_CI: is_trustworthy = True
        else:
            if dist < self.trust_radius_general: is_trustworthy = True
        #
        return request, is_trustworthy

    def loadweights(self, filename):
        """ Weights are loaded for the interpolators from a file. As the
            NearestNeighborInterpolator is not using the save option, also
            here, interpolators are just set up from the database

            Parameters:
            -----------
                filename, str:
                    filepath of the file containing the weights. Not used here!
        """
        #
        self.logger.warning("NearestNeighborInterpolator cannot load weights, interpolators are " +
                            "set up from DB")
        # As saving is not used, interpolators are set up from the database
        self.get_interpolators(self.db, self.properties)

    def save(self, filename):
        """ Method to save the interpolators to a file. Not used here!
            
            Parameters:
            -----------
                filename:
                    filepath where to save the information. Not used here!
        """
        #
        self.logger.warning("NearestNeighborInterpolator cannot be saved to a file")

#    @Timer(name="train")
    def _train(self):
        """ Method to train the interpolators. In the case of the NearestNeighborInterpolator
            the cKDTree has to be updated.
        """
        #update cKDTree
        self.tree = cKDTree(self.crds)

class NNInterpolator():
    """ NearestNeighborInterpolator for one property. """
    def __init__(self, db, ckdtree, prop, norm=2):
        """ 
            Parameters:
            -----------
                db: 
                    database containing the datasets on which the interpolation is based on

                ckdtree:
                    cKDTree of the coordinates for the interpolation

                prop: str
                    property that should be fitted. No sanity check with the database
                    is made!

                norm: 
                    norm for the cDKTree.query
        """
        #
        self.db = db
        self.tree = ckdtree
        self.prop = prop
        self.norm = norm

    def __call__(self, crd, request, idx=None):
        """ Returns the desired property of the nearest neighbor to the given geometry
            
            Parameters:
            -----------
                crd:
                    coordinates where the property is requested. The shape has to be consistent
                    with the shape of the coordinates in the cKDTree

                request:
                    Instance of the SPP request. Not used here, but needed for consistency.

                idx: int, optional
                    if idx is an integer, no nearest neighbor search is performed, but the 
                    property of the specific index is returned

            Returns:
            --------
                entry of the DB of the property next to the desired coordinate
        """
        #
        if idx is None:
            dist, idx = self.tree.query(crd, p=self.norm)
        return self.db.get(self.prop, idx)
