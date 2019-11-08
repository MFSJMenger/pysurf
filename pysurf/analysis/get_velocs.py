import sys
import os 
import numpy as np
import click

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.utils.chemutils import get_atom_from_mass

from pysurf.utils.constants import bohr2angstrom

def write_veloc(atoms, coord, step):
    string = str(len(coord)) + '\n'
    string += 'step {0} \n'.format(step)
    for i in range(len(coord)):
        string += '{0:s}  {1:12.8f}  {2:12.8f}  {3:12.8f}\n'.format(atoms[i], *coord[i])
    return string

def write_veloc_model(coord, step):
    string = str(len(coord)) + '\n'
    string += 'step {0} \n'.format(step)
    np.vectorize(str)
    string += np.array2string(coord, separator=',   ', precision=5).strip(']').strip('[')
    string += '\n'
    string += '\n'
    return string

@click.command()
@click.option('-f', 'infile', default='prop.db')
@click.option('-o', 'outfile', default='veloc.dat')
def get_velocs_command(infile, outfile):
    get_velocs(infile, outfile)

def get_velocs(infile, outfile):
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
        for m in range(len(db['coord'][0])):
            atoms+=['Q']
    if model is False:
        for m in mass[:,0]:
            atoms += [get_atom_from_mass(m)]
    
    
    with open(outfile, 'w') as output:
        step = 0
        for veloc in db['veloc']:
            if model is False:
                output.write(write_veloc(atoms, veloc, step))
            else:
                output.write(write_veloc_model(veloc, step))
            step += 1
        
if __name__=="__main__":
    get_velocs_command()
