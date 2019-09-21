from abc import ABC, abstractmethod


class SurfaceHoppingBase(ABC):


    def __init__(self, configfile, SPP):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model

        """
        pass

