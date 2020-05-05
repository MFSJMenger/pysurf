========
SAMPLING
========

Sampling:
---------

The sampling module allows to generate sampling for spectral or dynamical calculations. Different
samplers are added as Colt Plugins. There are two different types of samplers. Samplers only providing
coordinates (inherited from CrdSamplerBase) and samplers providing coordinates and velocities
(inherited from DynSamplerBase) as needed for dynamical simulations. Accordingly, there are two different
kinds of conditions, one with coordinates only (CrdCondition) and one with coordinates and velocities (DynCondition).





SamplingDB:
-----------

The sampling module comes with its own database representation. It is a PySurf Database, but a few more methods are added.

