import os
import numpy as np
#
from pysurf.colt import Colt
from pysurf.colt import FromCommandline
#
from pysurf.utils.constants import fs2au
from pysurf.sampling.initialconditions import InitialConditionsFactory


@FromCommandline("""
inputfile = initconds.inp :: file 
""")
def command_setup_initconds(inputfile):
    """ Setting up initial conditions according to the inputfile.
    If inputfile doesn't exist, colt will ask all necessary questions
    """
    initconds = InitialConditions.from_inputfile(inputfile)

if __name__=="__main__":
    command_setup_initconds()
