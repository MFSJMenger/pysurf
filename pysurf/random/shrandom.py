from abc import ABC, abstractmethod
import numpy as np


# This should be a Singleton class!
# there can only exist ONE way to 
# generate random numbers!

class RandomNumberGeneratorBase(ABC):
    """"Python Class to generate random numbers"""

    def __init__(self, seed=None):
        self.seed = seed
        self._set_seed(seed)

    @abstractmethod
    def _set_seed(self, seed):
        """Sets seed for random number generator"""
        pass

    @abstractmethod
    def get(self):
        pass


class RandomNumberGeneratorNP(RandomNumberGeneratorBase):

    def _set_seed(self, seed):
        """Sets seed for random number generator"""
        np.random.seed(seed)

    def get(self):
        return np.random.random()
