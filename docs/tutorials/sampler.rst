Tutorial: Write a Sampler
-------------------------

In case you want to extend the sampling facilities its rather
simple to write an own **Sampler**.

First thing you have to decide is what kind of sampler you want to
write. PySurf comes with two base *Sampling* classes:

    1. CrdSampler:
      Most basic sampler, it only returns coordinates of the system.
         - **TODO:** add features

   2. DynSampler:
      Sampler used for *dynamics* simulations, it needs to return
      **initial condition**, meaning positions and velocities.
         - **TODO:** add features
      All DynSampler are also at the same time CrdSampler.


Example: Gromacs Trajectory Sampler
-----------------------------------

In this tutorial, we show how to write a simple DynSampler,
that uses MDAnalysis_ to read information from a Gromacs Trajectory
and reads both coordinates and velocities.

Step 0: Basic Features
~~~~~~~~~~~~~~~~~~~~~~

As we are writing a *DynSampler*, we need to import **DynSamplerBase**.
Additionally, we will use numpy an MDAnalysis_ in this tutorial, so
we start with:

.. code-block:: python
   :linenos:

   import numpy as np
   import MDAnalysis as mda
   from pysurf.sampler import DynSamplerBase

Step 1: Basic Features
~~~~~~~~~~~~~~~~~~~~~~

Every **DynSamplerBase** has to implement 3 functions:

   from_config(cls, config), classmethod:
      classmethod to initialize the class using the information
      provided from the user

   get_condition(self):
      returns the next condition

   get_init(self):
      returns the initial structure


.. code-block:: python
   :linenos:

   class GromacsSampler(DynSamplerBase):
      """Sample from Gromacs Trajectory"""

      @classmethod
      def from_config(cls, config):
         """Initialize the object"""
         return cls()

      def get_condition(self):
         """return the next condition"""

      def get_init(self):
         """return the initial condition"""

Step 2: User Input
~~~~~~~~~~~~~~~~~~

PySurf is build around the Colt_ framework, developed along the lines
of this project. To specify specific input needed for your class you simply
use the *_questions* string:
In our example we will need 4 user inputs:

   - tpr_file, existing_file:
      the **tpr** file of Gromacs
   - trr_file, existing_file:
      the **trr** file of Gromacs
   - every, int:
      read only every nth step
   - selection, str, optional:
      save only a specific MDAnalysis-Selection_, if not set, use all atoms

.. code-block:: python
   :linenos:

   class GromacsSampler(DynSamplerBase):
      """Sample from Gromacs Trajectory"""

      _questions = """
      # name of the tpr file
      tpr_file = :: existing_file
      # name of the trr file
      trr_file = :: existing_file
      # read every nth trajectory from the file
      every = 1 :: int :: >1
      # selection
      selection = :: str, optional
      """

      @classmethod
      def from_config(cls, config):
         return cls(config['tpr_file'], config['trr_file'], every=config['every'],
                    selection=config['selection'])
 
      def __init__(self, tpr_file, trr_file, every=1, selection=None):
         # MDA Universe
         self.universe = mda.Universe(tpr_file, trr_file)
         # Maximum number of frames
         self.n_max = len(self.universe.trajectory)
         # current active frame
         self.current = 0
         # save every
         self.every = every
         # get the selection
         if selection is None:
            # use all atoms
            self.atoms = self.universe.atoms
         else:
            # use a specific selection
            self.atoms = self.universe.select_atoms(selection)


Step 3: Get Condition
~~~~~~~~~~~~~~~~~~~~~
The next and last step is to implement the get_condidion function, to
actually return the correct information back the the sampler:

.. code-block:: python
   :linenos:

   class GromacsSampler(DynSamplerBase):
      """Sample from Gromacs Trajectory"""

      ...

      def get_condition(self):
         if self.current < self.n_max:
            return self.next_condition()
         raise IndexError(f"Trajectory has only {self.n_max} elements, cannot access {self.current}")

We check that the current trajectory is available, and if so, we read it and move to the next
condition. If not we raise an **IndexError**.

.. code-block:: python
   :linenos:

   class GromacsSampler(DynSamplerBase):
      """Sample from Gromacs Trajectory"""

      ...

      def get_condition(self):
         if self.current < self.n_max:
            return self.next_condition()
         raise IndexError(f"Trajectory has only {self.n_max} elements, cannot access {self.current}")

      def next_condition(self):
         crd = np.copy(self.atoms.positions)
         veloc = np.copy(self.atoms.velocities)
         # 
         self.update_condition()
         return self.condition(crd=crd, veloc=veloc, state=None)

      def update_condition(self):
         try:
            for _ in range(self.every):
               next(self.universe.trajectory)
            self.current += self.every
         except StopIteration:
            self.current = self.n_max


.. _MDAnalysis: https://www.mdanalysis.org/
.. _MDAnalysis-Selection: https://www.mdanalysis.org/docs/documentation_pages/selections.html
.. _Colt: https://github.com/mfsjmenger/colt

