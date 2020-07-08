Tutorial: Trajectory Surface Hopping
====================================

In this tutorial we will show how to perform a trajectory surface hopping calcualation
within the `PySurf` framework. As an example we choose one of the models and put the 
trajectories on the second excited state. Finally we analyse the results.


Step 1: Sampling
================

First we have to generate the initial conditions for the propagation. PySurf comes
with several samplers and new ones can always be added. For a propagation it is important
to use a so-called `DynSampler`, which provides coordinates and velocities. For our
example we use the `Wigner` sampler. Make a folder where you want to perform the dynamics and
to do the initial sampling type:

.. code-block:: bash
    
    mkdir dynamics
    cd dynamics
    python <pysurf>/bin/sampling.py

<pysurf> stands for the path where `PySurf` is installed. Then we have to fill in the questions.
For educational purposes we only use 10 trajectories. Therefor answer:

.. code-block:: bash
    n_conditions = 10
    method = Wigner
    sampling_db = sampling.db
    model = True
    ..
    model = PyrazineSala
    n_states = 3

This part has to be rewritten with the new Sampler....

PySurf generates a database called sampling.db in the folder, that contains the conditions.


Step 2: Setup the calculations
==============================

Now we setup the calculations for each trajectory using the `setup_prop_calc.py` script.

.. code-block:: bash
    python <pysurf>/bin/setup_propagation.py 2

The positional argument of the script is the state, where the trajectories should start from.
Note that in PySurf the ground-state is 0, the first excited state is 1 and the second excited state is two 
and so on. The number of states includes the ground-state. In this case we have 3 states in total and start from 
the second excited state.

The script also collects the information for the SPP. All corresponding questions will be asked. In our case, we 
state that we want to get the information from a model, i.e. the PyrazineSala model with 3 states.

The `setup_propagation.py` script creates a folder named `prop` and in that folder there are folders for each trajectory
starting with `traj_` and then a 8 digit number.  


Step 3: Start the calculations
==============================

Pysurf offers comfortable tools to solve all tasks that have to do with the folder structure. In principle the
PySurf folder structure is that you always have a main folder and subfolders. The subfolders start with a common name,
e.g. `traj_` followed by a number. When you look into your `prop` folder, which has been created by the `setup_propagation.py`
script, there are 10 folders with `traj_00000000` up to `traj_00000009`. 

