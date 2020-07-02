Tutorial: Add a Model System
----------------------------

It is very easy and straightforward to add your own model system
to the PySurf Framework via the Plugin engine. In the models,
specific user input can be put in and your model is natively added
to the program package. Just put your python script in the plugin folder.


Example:  1D Harmonic Oscillator
--------------------------------

In this tutorial, we will implement a simple model with one or two harmonic 
oscillators.The user can decide, whether the second surface is included and whether
it should be shifted, as well as the frequencies and the offset energy.


Step 0: Colt Factory
~~~~~~~~~~~~~~~~~~~~~~

*Models* are written as Plugins of the *Model* 
Colt-Factory, i.e. they are inherited from the Factory. 


.. code-block:: python
   :linenos:

   from pysurf.spp import Model
   from pysurf.system import Mode


Step 1: Basic Features
~~~~~~~~~~~~~~~~~~~~~~

Every **Model** has to provide some features:

    implemented
        which properties are provided by the model, e.g. ['energy', 'gradient', 'fosc']

    modes
        list of the ground state modes for the Samplers

    @classmethod
    def from_config
        takes user input if questions are there and has to call the init function

    def get(self, request):
        fill the request with the demanded data 


.. code-block:: python
   :linenos:

    class HarmonicOscillator1D():
        """ Model for a 1D harmonic oscillator with 1 or 2 PES """

        implemented = ['energy', 'gradient']
        crd = [0.0]
        frequencies = [1.0]
        displacements = [[1.0]]
        modes = [Mode(freq, dis) for freq, dis in zip(frequencies, displacements)]

        @classmethod
        def from_config(cls, config):
            """initialize class with a given questions config!"""

        def get(self, request):
            """get requested info"""



Step 2: User Input
~~~~~~~~~~~~~~~~~~

PySurf is build around the Colt_ framework, developed along the lines
of this project. To specify specific input needed for your class you simply
use the *_questions* string:
In our example we will need user input for the shape of the potentials:
   - E0/E1
      offset energy of the corresponding potential
   - x0/x1
      shift of the harmonic oscillator in x direction 
   - w0/w1
      frequency of the corresponding PES
   - npes
      number of PES. The parameters for the second PES are only asked if
      npes = 2


.. code-block:: python
   :linenos:

    class HarmonicOscillator1D(Model):
    """ Model for a 1D harmonic oscillator with 1 or 2 potential energy surfaces """

    _questions = """
    e0 = 0.0 :: float
    w0 = 1.0 :: float
    x0 = 0.0 :: float
    # Number of potential energy surfaces
    npes = 1 :: str ::

    [npes(1)]

    [npes(2)]
    e1 = 1.0 :: float
    w1 = 1.0 :: float
    x1 = 1.0 :: float 
    """

    implemented = ["energy", "gradient"]
    masses = [1.0]
    crd = [0.0]
    frequencies = [1.0]
    displacements = [[1.0]]
    modes = [Mode(freq, dis) for freq, dis in zip(frequencies, displacements)]

    @classmethod
    def from_config(cls, config):
        e0 = config['e0']
        w0 = config['w0']
        x0 = config['x0']
        npes = int(config['npes'].value)

        config_npes = config['npes']
        return cls(e0, w0, x0, npes, config_npes)

    def __init__(self, e0, w0, x0, npes, config_npes):
        self.frequencies = [w0]
        self.crd = [x0]
        self.npes = int(npes)
        self.w = [w0]
        self.x = [x0]
        self.e = [e0]

        if self.npes == 2:
            self.w += [config_npes['w1']]
            self.x += [config_npes['x1']]
            self.e += [config_npes['e1']]





According to the Colt style, all arguments in the main question block are given in the
init method explicitely, whereas other question blocks are passed as config of the block.




Step 3: Implement *_energy* and *_gradient* function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The next step is to implement the actual model, i.e. the properties of the model.
In our case we want to provide the energies and the gradients of the surfaces and thus
private methods for *_energy* and *_gradient* are implemented.

.. code-block:: python
   :linenos:

    class HarmonicOscillator1D(Model):
    """ Model for a 1D harmonic oscillator with 1 or 2 potential energy surfaces """

      ...
      
        def _energy(self, x):
            energy = []
            for i in range(self.npes):
                energy += [0.5*self.w[i]*(x - self.x[i])**2 + self.e[i]]
            energy = np.array(energy).flatten()
            return energy

        def _gradient(self, x):
            gradient = {}
            for i in range(self.npes):
                gradient[i] = np.array(self.w[i]*(x - self.x[i]))
            return gradient


The *_energy* function takes a coordinate position, i.e. a numpy array and
returns an array with one or two entries, depending on whether the
surface contains one or two potential energy surfaces. The energy is calculated
according to the formula $0.5*\omega*(x-x0)^2$

The *_gradient* function takes a coordinate position, i.e. a numpy array and
returns a dictionary. The keys in the dictionary are the state numbers as integers,
i.e. 1 or 2 and the values are the gradients of the corresponding state.
In our specific case, such a dictionary may look like: {1: [0.5], 2: [0.2]}



Step 4: Implement the *get* function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The *get* function is called with the *request* as parameter. It 
has to fill in the desired results from the model into the *request* instance.

.. code-block:: python
   :linenos:

    class HarmonicOscillator1D(Model):
    """ Model for a 1D harmonic oscillator with 1 or 2 potential energy surfaces """

    ...

    def get(self, request):
        """the get function returns the adiabatic energies as well as the
           gradient at the given position crd. Additionally the masses
           of the normal modes are returned for the kinetic Hamiltonian.
        """
        crd = request.crd
        print('crd', crd)
        for prop in request:
            if prop == 'energy':
                request.set('energy', self._energy(crd))
            if prop == 'gradient':
                request.set('gradient', self._gradient(crd))
        return request



The *get* function checks which properties are demanded in the request and fills them into
the request, using the *_energy* and *_gradient* functions.





.. _Colt: https://github.com/mfsjmenger/colt

