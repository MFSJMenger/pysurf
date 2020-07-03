Tutorial: Write your own Workflow
=================================

Every scientist is after their own ideas and visions to move the boarders of knowledge.
Therefor, it is natural that predefined tools can never be sufficient for all tasks that scientists
want to do. The PySurf Workflow engine provides a toolbox of nodes, which can be combined like
"Lego Bricks" to new powerful algorithms. If your desired functionality is not yet in the toolbox,
you can easily add it and you can of course include the full functionality of the Workflow nodes
in your python scripts. The Workflow engine comes with its own domain specific language to combine 
all different nodes like in a normal script. The engine checks that the input and output of the nodes
which are combined fit together.
Let's see how it works in an example!

Example:  Single Point Calculation
----------------------------------

In this tutorial, we will implement a Workflow for a single point calculation and
print the results.


Step 0: Workflow framework
~~~~~~~~~~~~~~~~~~~~~~~~~~

For the Workflow framework, the *engine* has to be imported from *pysurf.workflow*.
With the command *engine.create_workflow()* the workflow is generated. The first argument
is a name, the second argument is a multiline string which contains the workflow.
With the command *workflow.run()* the workflow is executed.


.. code-block:: python
   :linenos:

    from pysurf.workflow import engine
    
    
    workflow = engine.create_workflow("populations", """
    ...
    """)
    
    wf = workflow.run()


Step 1: Include Workflow Nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. code-block:: python
   :linenos:
    
   workflow = engine.create_workflow("sp_calc", """
   crd = read_xyzfile_crd(crd_file)
   atomids = read_xyzfile_atomids(crd_file)
   spp = spp_calc("spp.inp", atomids, nstates, properties=properties)
   res = sp_calc(spp, crd, properties=properties)
   """)


- read_xyzfile_crd:
   node that returns the xyz coordinates from a xyz file
- read_xyzfile_atomids
  node that returns the atomids from a xyz file
- spp_calc
  node that initializes a SPP using an inputfile (filename, *file*), atomids (integer list, *ilist*), number of states (integer, *int*) and the desired
  properties (*list*)
- sp_calc
  node that sends the request to an initialized SPP. To start a calculation it is important that the SPP has been initialized with *spp_calc* and not
  *spp_analyse*. The second uses interpolation to produce the results. As arguments it takes the SPP (*spp*), the coordinates (*crd*) and the properties
  (*list*)

Variables which are not defined within the workflow are asked via the command line. In this case the user has to specify the xyz file (*crd_file*), the
number of states (*nstates*) and the properties that should be calculated as list, e.g. ['energy', 'gradient', 'fosc']

Step 2: Using the results
~~~~~~~~~~~~~~~~~~~~~~~~~

All variables of a workflow a saved in a dictionary. Either results are used within the workflow and saved to files etc., or 
you can read them from the dictionary and used them in your own script. Of course the variable *res* has to contain *energy*!

.. code-block:: python
   :linenos:
    
    workflow = engine.create_workflow("populations", """
    ...
    res = sp_calc(spp, crd, properties=properties)
    """)
    
    wf = workflow.run()
    print(wf['res']['energy'])



Appendix: How to put data into the Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is not only possible to take data out of the workflow, but also to put it in via a dictionary
when the workflow is executed.


.. code-block:: python
   :linenos:
    
    workflow = engine.create_workflow("populations", """
    ...
    """)

    wf = workflow.run({"properties": ['energy']})




.. _Colt: https://github.com/mfsjmenger/colt

