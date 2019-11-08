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

      "Request"  -> "SPP";
      "SPP" -> {"model", "DB", "Interface"};
   }

Hereby, it can get its Information potentially 3 ways
(Current implementation supports 0 xD), namely:

    1. from a defined model
    2. with the database 
    3. directly from an electronic structure

.. autoclass:: pysurf.wigner.WignerSampling
   :members:


Database
--------

.. graphviz::

   digraph database {

     node [shape=record];


      "DB"  -> "Interface" -> "Electronic Structure";

   }

Interface
---------

.. autoclass:: pysurf.wigner.WignerSampling
   :members:

