import os
import numpy as np
import click
from shutil import copy2

from pysurf.colt import Colt
from pysurf.sh.run_trajectory import RunTrajectory
from pysurf.utils.osutils import exists_and_isfile
from pysurf.logger import get_logger
from pysurf.initconds.initialconditions import InitialConditions




@click.command()
@click.option('-f', 'filename', default='propagation.inp')
def command_run_trajectory(filename):
    RunTrajectory(filename)

if __name__=="__main__":
    command_run_trajectory()
