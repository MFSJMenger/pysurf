Tutorial: Write an Interpolator
-------------------------------

In case you want to extend the interpolation facilities its
straightforward to write an additional **Interpolator**.

First, it has to be decided what interpolator should be
implemented. In this tutorial a nearest neighbor interpolator
is implemented.


Example:  Nearest Neighbor Interpolator
---------------------------------------

A Nearest Neighbor interpolator just looks for the dataset, whose
distance in the coordinate space is smallest. Then, the properties
of this nearest dataset are returned for the requested coordinates.
Our implementation in this tutorial will use scipy's nearest neighbor
search algorithm for a fast search and overall performance

Step 0: Colt Factory
~~~~~~~~~~~~~~~~~~~~~~

*Interpolator* are written as Plugins of the *Interpolator* 
Colt-Factory, i.e. they are inherited from the Factory. 
The Factory takes already care of many things, which are common
to any interpolator, e.g. coordinate transformation and 
energy-only calculations etc. During the example we will 
encounter many of them and explain them whenever they appear.


.. code-block:: python
   :linenos:

   from scipy.spatial import cKDTree 
   from pysurf.spp import Interpolator

Step 1: Basic Features
~~~~~~~~~~~~~~~~~~~~~~

Every **Interpolator** has to implement 5 functions:

    def get(self, request):
        fill the request with the demanded data and additionally
        a flag has to be returned whether the result is trustworthy
        If *fit_only* of the SPP is False, an electronic structure
        calculation is started, else the result of the interpolator 
        is taken. If *fit_only* is True, the result of the interpolator
        will be used anyway.

    def get_interpolators(self, db, properties):
        for each property that is requested, a separate interpolator is 
        set up. The function has to return a dictionary that contains the 
        property name and the interpolator for the property

    def save(self, filename):
        The save method is primarily for machine learning algorithms. The
        weights have to put in a file.

    def get_interpolators_from_file(self, filename, properties):
        This function reads the weightsfile and sets up the interpolators
        for the properties from the weightsfile. Specifically in 
        machine learning algorithms, its not necessary to train the 
        algorithm again.

    def _train(self):
        This function uses the input data to train the interpolator. It is
        primarily used in machine learning interpolators.

    def loadweights(self, filename):
        The weights are loaded from a file.



.. code-block:: python
   :linenos:

    class NearestNeighborInterpolator(Interpolator):
        """Nearest Neighbor Interpolator"""

        @classmethod
        def from_config(cls, config, db, properties, logger, energy_only, weightsfile, crdmode,
                        fit_only):
            """ This classmethod overwrites the from_config method of the factory and is needed
                if interpolator specific user input is implemented
            """

        def __init__(self, db, properties, logger, energy_only=False, weightsfile=None, 
                     crdmode='cartesian', fit_only=False)
            """ This init overwrites the init of the interpolator factory. Therefor, it is
                advisable to call the init of the interpolator factory to make sure that all
                functionality is working and only additional features are added here.
            """
            super().__init__(db, properties, logger, energy_only, weightsfile, crdmode=crdmode,
                             fit_only=fit_only)


        def get(self, request):
            """fill request

               Return request and if data is trustworthy or not
            """

        def get_interpolators(self, db, properties):
            """ """

        def save(self, filename):
            """Save weights"""

        def get_interpolators_from_file(self, filename, properties):
            """setup interpolators from file"""

        def _train(self):
            """train the interpolators using the existing data"""

        def loadweights(self, filename):
            """load weights from file"""


Step 2: User Input
~~~~~~~~~~~~~~~~~~

PySurf is build around the Colt_ framework, developed along the lines
of this project. To specify specific input needed for your class you simply
use the *_questions* string:
In our example we will need 4 user inputs:
   - trust_radius_general, float
      the radius to decide whether an interpolation is trustworthy
   - trust_radius_ci, float
      the radius in the region of small energy gaps to decide
      whether an interpolation is trustworthy
   - energy_threshold, float
      the threshold to distinguish between regions with small and
      large energy gap
   - norm, str 
      the norm that is used to measure the distance between points


.. code-block:: python
   :linenos:

    class NearestNeighborInterpolator(Interpolator):
        """Basic Rbf interpolator"""
    
        _questions = """
            trust_radius_general = 0.75 :: float
            trust_radius_ci = 0.25 :: float
            energy_threshold = 0.02 :: float
            norm = euclidean :: str :: [euclidean]
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
            self.trust_radius_general = trust_radius_general
            self.trust_radius_CI = trust_radius_CI
            self.energy_threshold = energy_threshold
            self.tree = None
            self.norm = norm
            # Call the init method of the Interpolator Factory
            super().__init__(db, properties, logger, energy_only, weightsfile,
                             crdmode=crdmode, fit_only=fit_only)


Parameters

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



Step 3: Implement *get_interpolators* function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The next step is to implement the *get_interpolators* method and the helper class
for the NearestNeighborInterpolator of each property *NNInterpolator*. For each
property, an Interpolator is set up, which is an instance of the *NNInterpolator*
class. Each interpolator has to be callable and to return the desired property.

.. code-block:: python
   :linenos:

   class NearestNeighborInterpolator(Interpolator):
      """Nearest Neighbor Interpolator"""

      ...
        def get_interpolators(self, db, properties):
            """ """
            self.tree = cKDTree(self.crds)
            return {prop_name: NNInterpolator(db, self.tree, prop_name)
                    for prop_name in properties}, len(db)
    

    class NNInterpolator():
        def __init__(self, db, ckdtree, prop):
            self.db = db
            self.tree = ckdtree
            self.prop = prop

        def __call__(self, crd, request=None, idx=None):
            if idx is None:
                dist, idx = self.tree.query(crd)
            return self.db.get(self.prop, idx)

The *get_interpolators* method returns a dictionary with the property names as keys and the
interpolator for that specific property as value. For each property a separate interpolator
has to be set up so that the interpolator factory can handle the interpolators for the properties
independently, which allows e.g. the energy_only calculations. Implementing the interpolators in
this way, they naturally are included in the code package and the full functionality is available.


To avoid that the cKDTree is set up several times, the *NNInterpolator* takes the
tree as a Parameter. Moreover, if *NNInterpolator* is called with an index,
no nearest neighbor search is performed, but the property of dataset with the index
is returned. This is important in the case when several properties are demanded so that the
nearest neighbor search is done only once, cf. Step 4 and the *get* function.




Step 4: Implement *get* function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The *get* function is called with the *request* as parameter. It 
has to fill in the desired results from the fit into the *request* instance and state
whether the fit is trustworthy.

.. code-block:: python
   :linenos:

   class NearestNeighborInterpolator(Interpolator):
      """Nearest Neighbor Interpolator"""

      ...

    def get(self, request):
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


The *get* function first has to make sure that the interpolators get the right coordinates.
Subsequently, the interpolators for all the desired properties are called and the
results are put into the request instance. 
Finally, it is checked, whether the requested point is within the trusted region. The trusted
region is devided into two parts, depending whether the smallest energy gap between two potential
energy surfaces is small or large. The threshold is given as the *energy_threshold* as user input
as well as the radii *trust_radius_ci* and *trust_radius_general*.



Step 5: Implement the save, load and _train methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The NearestNeighborInterpolator does not use the functionality that the interpolators and their parameters
are stored to a file and read afterwards, to avoid long training sessions. The training of the Nearest Neighbor
search is just to update the cKDTree, which doesn't take very long. Therefor, these functions are not really used,
but implemented in a way to make sure that the full functionality is available.



.. code-block:: python
   :linenos:

   class NearestNeighborInterpolator(Interpolator):
      """Nearest Neighbor Interpolator"""

      ...

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

    def _train(self):
        """ Method to train the interpolators. In the case of the NearestNeighborInterpolator
            the cKDTree has to be updated.
        """
        #update cKDTree
        self.tree = cKDTree(self.crds)



.. _Colt: https://github.com/mfsjmenger/colt

