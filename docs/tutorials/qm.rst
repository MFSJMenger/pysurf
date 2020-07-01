Tutorial: Write a Abinitio-Interface
====================================

Abinitio-Interfaces are `plugins` for the `SurfacePointProvider` and are used 
to perform (electronic structure) calculations on molecular systems.
In the following example implementations are presented in order to 

Example: Adopter to the Atomic Simulation Environment
-----------------------------------------------------
According to the documentation, the Atomic Simulation Environment (ASE_) is `a set of tools and
Python modules for setting up, manipulating, running, visualizing and analyzing atomistic
simulations` under GNU LPGL licence. 
Hereby, the ASE provides so called `Calculators`, which are similar to the Abinitio-Interfaces used in 
the `SurfacePointProvider`.
In the following, we are going to construct a adopter class to use the ASE calculators within out
framework.

Step 0: Basic Setup
~~~~~~~~~~~~~~~~~~~
To implement a Abinitio-Interface using the ASE we need to import some
objects/methods.
We will use from the ASE: 

    - `Atoms`: which defines the basic molecule
    -`calculators`: which gives access to the available calculators

From PySurf we will import the Abinitio class which is the baseclass 
of all Abinitio interfaces


.. code-block:: python
   :linenos:

   from ase import Atoms
   from ase import calculators
   from pysurf import Abinitio


Step 1: Abstractmethods in Abinitio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. code-block:: python
   :linenos:

    class ASEInterface(Abinitio):
        implemented = ['energy', 'gradient']

        @classmethod
        def from_config(cls, config, atomids, nstates):
            """used to initialize the class using userinput"""

        def get(self, request):
            """Compute requested properties and set results"""

The `Abinitio` baseclass defines two abstractmethods that need to be set:

    - from_config: 
        used to initialize the class using the userinput
    - get:
        which is used to answer the request
    implement:
        states which properties are implemented, for runtime check
        do not lie on that...

If you write a new qm interface, those are the methods you have to implement.

The next thing is to add user configurations.


Step 2: Adding userinput for the Plugin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before we are going to implement the abstractmethods we are going to
add some custom user-input that we want to use in our `Plugin`. 
Herefor, we use the `question` DSL of Colt and add our own 
questions to our class.

.. code-block:: python
   :linenos:

    class ASEInterface(Abinitio):

        _questions =  """
        calculator = qchem 

        [calculator(qchem)]
        method = b3lyp
        basis = 6-31g

        [calculator(psi4)]
        method = b3lyp
        memory = 500MB
        basis = 6-31g
        """

        implemented = ['energy', 'gradient']

        @classmethod
        def from_config(cls, config, atomids, nstates):
            """used to initialize the class using userinput"""
            if nstates != 1:
                raise Exception("ASE does not support excited states")
            return cls(config['calculator'], atomids)

        def __init__(self, calculator, atomids):
            # define the molecule in a basic manner
            self.molecule = Atoms(numbers=atomids)
            self.calculator = self._select_calculator(calculator)

        def _select_calculator(self, calculator):
            if calculator == 'psi4':
                return calculators.psi4.Psi4(atoms=self.molecule, method=calculator['method'], 
                                             memory=Calculator['memory'], basis=calculator['basis'])
            if calculator == 'qchem':
                return calculators.qchem.QChem(atoms=self.molecule, method=calculator['method'], 
                                               basis=calculator['basis'])
            raise NotImplementedError("calculator not implemented")


For education purpose we only show two calculators, the one for qchem and the one for psi4.

Step 3: Implementing get
~~~~~~~~~~~~~~~~~~~~~~~~

With that it is now trivial to implement the get function

.. code-block:: python
   :linenos:

    class ASEInterface(Abinitio):
        ...

        def get(self, request):
            """Compute requested properties and set results"""
            # set the coordinates
            self.molecule.positions = request.crd
            # compute energy
            if 'energy' in request:
                request['energy'].set(self.calculator.get_forces())
            # compute gradient
            if 'gradient'in request:
                request['gradient'].set(self.calculator.get_forces())
            return request

With that done, we have a new Abinitio-Interface




.. _ASE: https://wiki.fysik.dtu.dk/ase/

