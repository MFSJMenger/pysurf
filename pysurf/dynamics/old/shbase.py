from abc import ABC, abstractmethod
from ..random.shrandom import RandomNumberGeneratorNP


class SurfaceHoppingBase(ABC):

    def __init__(self, configfile):
        """Setup surface hopping using config in `configfile`
        and a SurfacePointProvider (SPP) abstract class

        The system is completely defined via the SPP model
        """
        # get spp input
        spp_inp = self._parse_config(configfile) 
        # setup SPP
        self.spp = SurfacePointProvider(spp_inp)
        # setup trajectory
        self.irestart = self._init_trajectory()
        if self.irestart is True:
            self._restart_trajectory()
         
    @abstractmethod
    def _parse_config(self, config):
        """parse the config file and do anything related
           should return the input file for the spp
        """
        pass

    @abstractmethod
    def _restart_trajectory(self):
        """method to restart a trajectory"""
        pass

    @abstractmethod
    def _init_trajectory(self):
        """setup a new trajectory run, returns if restart or not!"""
        pass
    
    @abstractmethod
    def _init_system(self):
        """method to setup the surface hopping run"""
        pass

    @abstractmethod
    def run(self):
        """run the actual trajectory job!"""
        pass

    def init_random(self, seed=16661):
        self._random = RandomNumberGeneratorNP(seed)

    def random(self):
        return self._random.get()

    @staticmethod
    def x_update(crd, v, a, dt):
        """Currently we only support velocity verlet"""
        crd += v*dt + 0.5*a*dt*dt
        return crd

    @staticmethod
    def v_update(v, a_old, a, dt):
        v += 0.5 * (a_old + a) * dt
        return v

    @staticmethod
    def get_acceleration(g, m):
        return -g/m

    @staticmethod
    def _rescale_velocity_along_v(ekin, dE, v):
        factor = (1. - dE/ekin)
        v *= np.sqrt(factor)
        return v
