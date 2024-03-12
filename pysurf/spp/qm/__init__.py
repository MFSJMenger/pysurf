from ..methodbase import AbinitioBase
from .xtb import XTBInterface
from .qchem import QChem
from .turbomole import Turbomole
from .open_molcas import OpenMolcas 
from .bagel import Bagel 

plugins = [XTBInterface, QChem, Turbomole, OpenMolcas, Bagel]
