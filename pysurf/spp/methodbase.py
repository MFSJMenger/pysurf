from abc import abstractmethod
#
from .spp import AbinitioFactory, ModelFactory
# fileparser
from ..fileparser import read_geom


class AbinitioBase(AbinitioFactory):
    """Base Class for all Abinitio Implementations"""

    _register_plugin = False
    implemented = []

    @classmethod
    @abstractmethod
    def from_config(cls, config):
        """initialize class with a given questions config!"""
        pass

    @abstractmethod
    def get(self, request):
        """get requested info"""

    def _get_refgeo(self, filename):
        natoms, atoms, crds = read_geom(refgeo_path)
        return natoms, {'atoms': atoms, 'crd': np.array(crds)}

    def get_masses(self):
        # masses are given in an array of shape (natoms, 3) like
        # the coordinates so that they can be easily used in the
        # surface hopping algorithm
        return np.array([[[atomic_masses[atomname_to_id[self.refgeo['atoms'][i]]],
                    atomic_masses[atomname_to_id[self.refgeo['atoms'][i]]],
                    atomic_masses[atomname_to_id[self.refgeo['atoms'][i]]]]]
                    for i in range(self.natoms)
                ])


class ModelBase(ModelFactory):
    """Base Class for all Abinitio Implementations"""

    _register_plugin = False

    implemented = []

    @classmethod
    @abstractmethod
    def from_config(cls, config):
        """initialize class with a given questions config!"""
        pass

    @abstractmethod
    def get(self, request):
        """get requested info"""
        pass
