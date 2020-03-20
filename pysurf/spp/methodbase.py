from abc import abstractmethod
#
from .spp import AbinitioFactory, ModelFactory


class AbinitioBase(AbinitioFactory):
    """Base Class for all Abinitio Implementations"""

    _register_plugin = False

    @classmethod
    @abstractmethod
    def from_config(cls, config):
        """initialize class with a given questions config!"""
        pass

    @abstractmethod
    def get(self, request):
        """get requested info"""
        pass


class ModelBase(ModelFactory):
    """Base Class for all Abinitio Implementations"""

    _register_plugin = False

    @classmethod
    @abstractmethod
    def from_config(cls, config):
        """initialize class with a given questions config!"""
        pass

    @abstractmethod
    def get(self, request):
        """get requested info"""
        pass
