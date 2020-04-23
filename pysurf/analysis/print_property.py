#! /data/ehrmaier/anaconda3/bin/python3
import sys
import os 
import numpy as np

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.molecule.atominfo import get_atom_from_mass
from pysurf.colt import Colt

from pysurf.utils.constants import bohr2angstrom

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
        db = Database.load_db(config['infile'])
        
        for idx, entry in enumerate(db[config['property']]):
            print('\nEntry: {}'.format(idx))
            print(entry)
    
    @classmethod
    def from_config(cls, config):
        return cls(config)

if __name__=='__main__':
    PrintProperty.from_commandline()
