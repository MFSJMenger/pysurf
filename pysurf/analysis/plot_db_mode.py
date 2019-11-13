import numpy as np
import click

from pysurf.database.database import Database
from pysurf.molden import MoldenParser
from pysurf.wigner import WignerSampling
from pysurf.spp.spp import SurfacePointProvider

@click.command()
@click.argument('molden')
@click.argument('spp')
@click.argument('mode')
def plot_db_mode_command(molden, spp, mode):
    plot_db_mode(molden, spp, mode)

def plot_db_mode(molden, spp, mode):
    mode = int(mode)
    smp = WignerSampling.from_molden(molden)
    
    spp = SurfacePointProvider(spp)

    x = np.linspace(-5, 5, 100)
    y = np.zeros((100,spp.nstates), dtype=float)

    print('x', x)
    print('y', y)
    print('crd: ', smp.molecule.crd)
    print('mode: ', smp.modes[1])

    for c in x:
        crd = smp.molecule.crd + c * smp.modes[mode].displacements
        y = spp.get({'crd': crd, 'energy': None})

    print(y)


if __name__=="__main__":
    plot_db_mode_command()
