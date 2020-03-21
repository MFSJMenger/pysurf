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
from pysurf.colt import FromCommandline

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


@FromCommandline("""
infile = prop.db :: file_exists
outfile = crd.xyz :: file
""")
def get_coordinates_command(infile, outfile):
    get_coordinates(infile, outfile)

def get_coordinates(infile, outfile):
    if not(os.path.isfile(infile)):
        print('Error: infile path does not exist! ' + infile)
        exit()
    
    db = Database.load_db(infile)
    try:
        mass = db['mass']
        if len(mass.shape) == 1:
            model = True
        else:
            model = False
    except:
        print('Masses could not be found in DB!')
        mass = None
        model = True
    
    atoms=[]
    if model is True:
        for m in range(len(db['crd'][0])):
            atoms+=['Q']
    if model is False:
        for m in mass[:,0]:
            atoms += [get_atom_from_mass(m)]
    
    
    with open(outfile, 'w') as output:
        step = 0
        for crd in db['crd']:
            if model is False:
                output.write(write_xyz(atoms, crd, step))
            else:
                output.write(write_crd(crd, step))
            step += 1
    
if __name__=='__main__':
    get_coordinates_command()
