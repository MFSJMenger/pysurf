import sys
import os 
import numpy as np

from pysurf.database import PySurfDB
from pysurf.system.atominfo import get_atom_from_mass

from pysurf.utils.constants import bohr2angstrom
#
from colt import Colt

def au2ang(x):
    return x * bohr2angstrom

def write_xyz(atoms, crd, step):
    string = str(len(crd)) + '\n'
    string += 'step {0} \n'.format(step)
    for i in range(len(crd)):
        string += '{0:s}  {1:12.8f}  {2:12.8f}  {3:12.8f}\n'.format(atoms[i], *au2ang(crd[i]))
    return string

def write_crd(crd, step):
    string = str(len(crd)) + '\n'
    string += 'step {0} \n'.format(step)
    np.vectorize(str)
    string += np.array2string(crd, separator=',   ', precision=5).strip(']').strip('[')
    string += '\n'
    string += '\n'
    return string


class PrintProperty(Colt):

    _questions="""
        infile = :: existing_file

        property = :: str
        """

    def __init__(self, config):
        db = PySurfDB.load_database(config['infile'], read_only=True)

        if config['property'] == 'len':
            print("Number of entries in database: ", len(db))
            return

        for idx, entry in enumerate(db[config['property']]):
            print('\nEntry: {}'.format(idx))
            print(entry)
    
    @classmethod
    def from_config(cls, config):
        return cls(config)

if __name__=='__main__':
    PrintProperty.from_commandline()
