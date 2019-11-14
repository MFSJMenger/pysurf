import os
import numpy as np
import click

from pysurf.colt import AskQuestions
from pysurf.utils.constants import fs2au
from pysurf.initconds.initialconditions import InitialConditions

def initconds(inputfile):
    initconds = InitialConditions(inputfile)


@click.command()
@click.option('-f', 'inputfile', default='initconds.inp')
def command_initconds(inputfile):
    initconds(inputfile)

if __name__=="__main__":
    command_initconds()
