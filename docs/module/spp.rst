================
Surface Provider
================


SurfacePointProvider
--------------------

The SurfacePointProvider is an interface layer to get
the properties of a particular surface point defined
by the coordinates of a system.


.. graphviz::

   digraph surfacepointprovider {

     node [shape=record];

      "Request"  -> "SPP" -> "DB";
      "DB" -> {"Model", "Interface"};
      "SPP" -> {"Model", "Interface"};
   }

Hereby, it can get its information potentially 3 ways, namely:

    1. from a defined model
    2. directly from an electronic structure
    3. via interpolation from the database 


The main idea is the following: The SPP acts as a factory that
will create an general interface to any of the requested options.
It also checks that all possible requested features can be obtained
from the actual implementation. Meaning, if Gradients are needed, but
the underlying electronic structure method does not support that, then
the SPP will raise an exception.

Therefore the contract is the following, the caller asks the SPP what
setup to chose and then checks if all the possible request can be full-filled
by the chosen method. The SPP does not know directly any of the methods, 
but it knows how to setup the underlying interface and to get from there
all information it needs to do its check.

.. graphviz::

   digraph dependencies {

     node [shape=record];

      "SH"  -> "SPP" -> "Method";
   }


The SH Method is responsible to store information along the trajectory.
Hereby, the logging is partitioned into two kinds of data:

    1. Needed data in the run
    2. Additional, information (e.g. orbital coefficients, theodore analysis...)

Needed data, is logged according to the SH scheme, while all additional information
is directly logged after the SPP call.
The database needs to contain ANYTHING, that should be logged!


Database
--------

If a data

.. graphviz::

   digraph database {

     node [shape=record];


      "DB"  -> "Interface" -> "Electronic Structure";

   }

