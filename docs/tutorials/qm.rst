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




Example: PySCF interface:
-------------------------

In this example we show how to write an interface for the PySCF program 
package, which also supports excited states.

Step 0: Basic Setup
~~~~~~~~~~~~~~~~~~~
To implement a Abinitio-Interface using PySCF we need to import some
objects/methods.
We will use from the PySCF

    - `gto`: which defines the basic molecule
    - `dft`, `grad`, `tddft`: which allow to perform dft and tddft calculations for energies
      and gradients

From PySurf we will import the Abinitio class which is the baseclass 
of all Abinitio interfaces


.. code-block:: python
   :linenos:

   from pysurf import Abinitio
   from pyscf import gto, dft, tddft, grad


Step 1: Abstractmethods in Abinitio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. code-block:: python
   :linenos:

    class PySCF(Abinitio):

        methods = {}
        implemented = []

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

For each method a separate class is implemented, the calculator classes.
These calculator classes have to have methods with the name `do_prop` where
prop stands for all the implemented properties, e.g. `do_energy`. Moreover
it has to have a property `implemented` which is copied to the PySCF class.
PySurf will check the `implemented` property, whether the interaface provides
all necessary properties that are needed in the calculation.

.. code-block:: python
   :linenos:

    class PySCF(Abinitio):

        _questions =  """
        basis = 631g*
        # Calculation Method
        method = DFT/TDDFT :: str :: [DFT/TDDFT]
        """

        # implemented has to be overwritten by the individual classes for the methods
        implemented = []

        # dictionary containing the keywords for the method and the corresponding classes
        methods = {'DFT/TDDFT': DFT}

        @classmethod
        def _extend_questions(cls, questions):
            questions.generate_cases("method", {name: method.questions
                                                 for name, method in cls.methods.items()})

        @classmethod
        def from_config(cls, config, atomids, nstates):
            method = config['method'].value
            basis = config['basis']
            config_method = config['method']
            return cls(basis, method, atomids, nstates, config_method)


        def __init__(self, basis, method, atomids, nstates, config_method):
            """ """
            self.mol = self._generate_pyscf_mol(basis, atomids)
            self.nstates = nstates
            self.atomids = atomids
            self.basis = basis
            # initializing the class for the corresponding method
            self.calculator = self.methods[method].from_config(config_method, self.mol, nstates)
            # update the implemented property
            self.implemented = self.calculator.implemented


The code for the `_generate_pyscf_mol` function is shown in the next section. It is a PySCF specific
function that creates the molecule object for PySCF.

Step 3: Implementing get
~~~~~~~~~~~~~~~~~~~~~~~~

The get function calls the corresponding functions of the calculator class.
The `_generate_pyscf_mol` function generates the basic molecule object of Pyscf.

.. code-block:: python
   :linenos:

    class PySCF(Abinitio):
        ...
       # update coordinates
        self.mol = self._generate_pyscf_mol(self.basis, self.atomids, request.crd)
        for prop in request:
            func = getattr(self.calculator, 'do_' + prop)
            func(request, self.mol)
        #
        return request

        @staticmethod
        def _generate_pyscf_mol(basis, atomids, crds=None):
            """ helper function to generate the mol object for Pyscf """
            if crds is None:
                crds = np.zeros((len(atomids), 3))
            mol = gto.M(atom=[[atom, crd] for atom, crd in zip(atomids, crds)],
            basis = basis, unit='Bohr')
            return mol
    

Step 4: Implementing the DFT calculator class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For educational purposes we restrict to the calculation of energies.
Like in all `Colt` classes questions can be added, which are asked through
the `_extend_questions` method of the `PySCF` class. Answers are passed to 
the `__init__` function via the `from_config` classmethod.
At the initialization the dft and tddft scanners are set up to make sure
that calculations are started from the last converged result.
In the `do_energy` function the request is filled with the energies.

.. code-block:: python
   class DFT(Colt):
       """ class which executes the DFT and TDDFT calculations using the PySCF package """
   
       _questions = """
       functional = :: str :: ['pbe0']
       basis = ccpvdz :: str
       """
   
       implemented = ['energy']
   
   
       @classmethod
       def from_config(cls, config, mol, nstates):
           """ """
           functional = config['functional']
           return cls(functional, mol, nstates)
   
   
       def __init__(self, functional, mol, nstates):
           self.mol = mol
           self.nstates = nstates
   
           mydft = dft.RKS(mol).x2c().set(xc=functional)
           self.dft_scanner = mydft.as_scanner()

           if self.nstates > 1:
               # Switch to xcfun because 3rd order GGA functional derivative is not
               # available in libxc
               mydft._numint.libxc = dft.xcfun
               mytddft = tddft.TDDFT(mydft)
               self.tddft_scanner = mytddft.as_scanner()
               self.tddft_scanner.nstates = self.nstates - 1
    
    
       def do_energy(self, request, mol):
           if self.nstates == 1:
               en = [self.dft_scanner(mol)]
           else:
               en = self.tddft_scanner(mol)
    
           request.set('energy', en)



.. _PYSCF: https://sunqm.github.io/pyscf/

