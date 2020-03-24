import os
import numpy as np
import click

from pysurf.colt import Colt
from pysurf.colt import FromCommandline
from pysurf.utils.constants import fs2au
from pysurf.sampling.initialconditions import InitialConditionsFactory



#@FromCommandline("""
#inputfile = initconds.inp :: file 
#""")
@click.command()
@click.option('-f', 'inputfile', default='initconds.inp')
def command_setup_initconds(inputfile):
    initconds = InitialConditionsFactory.from_inputfile(inputfile)
#    for cond in initconds:
#        print(cond.veloc)

if __name__=="__main__":
    command_setup_initconds()
    #initconds.from_commandline(description="Basic Command line tool")
