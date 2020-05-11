from ..methodbase import AbinitioBase
from .xtb import XTBInterface
from .qchem import QChem
from .turbomole import Turbomole

plugins = [XTBInterface, QChem, Turbomole]
