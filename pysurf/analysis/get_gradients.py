import sys
import os 
import numpy as np

from pysurf.database.database import Database
from pysurf.database.dbtools import DatabaseRepresentation
from pysurf.database.dbtools import DatabaseTools
from pysurf.database.dbtools import DBVariable
from pysurf.database.dbtools import load_database
from pysurf.system.atominfo import get_atom_from_mass

from pysurf.utils.constants import bohr2angstrom
#
from colt import from_commandline


def write_gradient(atoms, gradient, step):
    nstates = len(gradient)
    string = str(len(gradient)) + '\n'
    string += 'step {0} \n'.format(step)
    for state in range(len(gradient)):
        string += 'state: ' + str(state) + '\n'
        for i in range(len(atoms)):
            string += '{0:s}  {1:12.8f}  {2:12.8f}  {3:12.8f}\n'.format(atoms[i], *gradient[state][i])

    return string

def write_gradient_model(gradient, step):
    string = str(len(gradient)) + '\n'
    string += 'step {0} \n'.format(step)
    np.vectorize(str)
    string += np.array2string(gradient, separator=',   ', precision=5).strip(']').strip('[')
    string += '\n'
    string += '\n'
    return string

@from_commandline("""
infile = prop.db :: file_exists
outfile = gradient.dat :: file
""")
def get_gradients_command(infile, outfile):
    get_gradient(infile, outfile)

def get_gradients(infile, outfile):
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
        for m in range(len(db['gradient'][0])):
            atoms+=['Q']
    if model is False:
        for m in mass[:,0]:
            atoms += [get_atom_from_mass(m)]
    
    
    with open(outfile, 'w') as output:
        step = 0
        for grad in db['gradient']:
            if model is False:
                output.write(write_gradient(atoms, grad, step))
            else:
                output.write(write_gradient_model(grad, step))
            step += 1
        
  
if __name__=="__main__":
    get_gradients_command()
