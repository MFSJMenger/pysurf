import sys
import os 
import numpy as np

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.utils.chemutils import get_atom_from_mass

from pysurf.utils.constants import au2ev



def write_energy(energy, step):
    string = ''
    string += 'step {0} :'.format(step)
    string += np.array2string(energy, separator=',   ', precision=5).strip(']').strip('[')
    string += '\n'
    return string


if len(sys.argv) < 2:
    print('Error: Please provide DB file')

infile = sys.argv[1]

if len(sys.argv) >=3:
    outfile = sys.argv[2]
else:
    outfile = 'energy.dat'

if not(os.path.isfile(infile)):
    print('Error: infile path does not exist! ' + infile)
    exit()

db = Database.load_db(infile)


with open(outfile, 'w') as output:
    step = 0
    for energy in db['energy']:
        output.write(write_energy(energy, step))
        step += 1
    

