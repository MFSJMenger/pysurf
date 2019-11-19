import os
import numpy as np
import click

from pysurf.colt import Colt
from pysurf.utils.constants import fs2au
from pysurf.sampling.initialconditions import InitialConditions

class initconds(Colt):
    _questions="""
    # Inputfile for initial conditions
    inputfile =
    """
    
    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, answer):
        initconds = InitialConditions.from_inputfile(answer['inputfile'])


if __name__=="__main__":
    initconds.from_commandline(description="Basic Command line tool")
