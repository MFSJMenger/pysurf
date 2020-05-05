import sys
import os 
import numpy as np

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.system.atominfo import get_atom_from_mass

from pysurf.utils.constants import au2ev
from pysurf.colt import FromCommandline



def write_etot(ekin, epot, etot, step):
    string = ''
    string += 'step {0} : {1:8.4f}, {2:8.4f}, {3:8.4f}'.format(step, ekin, epot, etot)
    string += '\n'
    return string


@FromCommandline("""
infile = prop.db :: file_exists
outfile = etot.dat :: file
""")
def get_etot_command(infile, outfile):
    get_etot(infile, outfile)

def get_etot(infile, outfile):
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
        
if __name__=="__main__":
    get_etot_command()
