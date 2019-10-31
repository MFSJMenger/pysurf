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



def write_etot(ekin, epot, etot, step):
    string = ''
    string += 'step {0} : {1:8.4f}, {2:8.4f}, {3:8.4f}'.format(step, ekin, epot, etot)
    string += '\n'
    return string


if len(sys.argv) < 2:
    print('Error: Please provide DB file')

infile = sys.argv[1]

if len(sys.argv) >=3:
    outfile = sys.argv[2]
else:
    outfile = 'etot.dat'

if not(os.path.isfile(infile)):
    print('Error: infile path does not exist! ' + infile)
    exit()

db = Database.load_db(infile)


with open(outfile, 'w') as output:
    output.write('ekin, epot, etot \n')
    step = 0
    for ekin, epot, etot in zip(db['ekin'], db['epot'], db['etot']):
        output.write(write_etot(ekin[0], epot[0], etot[0], step))
        step += 1
    

