Tutorial: Single Point Calculation
==================================

In this tutorial we will show how to perform an ab-initio single point calculation
using the `PySurf` framework. If you have any other electronic structure program installed
that is supported by PySurf, you can just use it. Otherwise you can download
the `PySCF` program package and install it via:

.. code-block:: bash

    pip install pyscf

`PySCF` is free of charge. In this tutorial we use it to do a TDDFT calculation.


Next you need a folder for the calculation and an xyz file. There is an example (where?) of a $SO_2$ geometry.
Go to your favorite place and make a folder for the test calculation and copy your coordinate file:

.. code-block:: bash

    mkdir test_so2
    cd test_so2
    cp <so2.xyz> ./

<so2.xyz> stands for the path of your so2 coordinate file.

Finally you have to run the `sp_calc.py` workflow for a single point calculation. <pysurf> is the path
where your PySurf package is installed

.. code-block:: bash

    python <pysurf>/bin/sp_calc.py so2.xyz 2 energy 

The `sp_calc.py` script takes three positional arguments. The first is the path of the coordinate file, the second is the number of states and the third is a
list of the properties that you want to calculate. If you also want to have the gradients you have to write:

.. code-block:: bash

    python <pysurf>/bin/sp_calc.py so2.xyz 2 "energy, gradient"

For help text of the script type:

.. code-block:: bash

    python <pysurf>/bin/sp_calc.py -h

Subsequently PySurf will guide you through all the questions of the `SurfacePointProvider` and the electronic structure
interface. You can take the defaults or choose your own. The `Colt` framework helps you that your answers are consistent.
In this case, we accept all the defaults up to the point for the `software`. There we choose `PySCF`. Subsequently we take 
the defaults again. The `PySCF` calculation is started and the output is printed to the screen. The results are saved in
the database that was put in the input (the default is just `db.dat`). To see all the parameters of your calculation,
you can open the `spp.inp` file, where the whole input is stored. If you want to repeat the calculation with the
same settings for a different geometry, you can just exchange the coordinate file. Typing again:

.. code-block:: bash

    python <pysurf>/bin/sp_calc.py so2.xyz 2 "energy, gradient"

Will start a new calculation. This time no questions will be asked, but all the answers are read from the `spp.inp` file.
See the other tutorials how to setup multiple calculations.

Congratualtions, you just performed your first single point calculation with PySurf!
