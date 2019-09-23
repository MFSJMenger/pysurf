from abc import ABC, abstractmethod
from ..random import RandomNumberGenerator


class SurfaceHoppingBase(ABC):


    def __init__(self, configfile, SPP):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model

        """
        pass

    @abstractmethod
    def _restart_trajectory(self):
        """method to restart a trajectory"""
        pass

    @abstractmethod
    def _init_trajectory(self):
        """method to initialize a trajectory"""
        pass
    
    @abstractmethod
    def init_system(self):

    def run_vv(iactive, nsteps, random_seed):

        random = RandomNumberGenerator(random_seed)

        for istep in range(nsteps):

