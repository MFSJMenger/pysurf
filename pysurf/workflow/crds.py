import numpy as np

from pysurf.spp import SurfacePointProvider
from pysurf.database import PySurfDB
from pysurf.system import ATOMNAME_TO_ID
from pysurf.qctools.converter import length_converter

from . import engine

@engine.register_action
def read_xyzfile(filename: "existing_file") -> "atomids_crd":
    conv = length_converter.get_converter('ang', 'au') 
    with open(filename, 'r') as infile:
        lines = infile.readlines()

    natoms = int(lines[0])
    atoms = []
    crds = []
    for i in range(natoms):
        atoms += [ATOMNAME_TO_ID[lines[2+i].split()[0]]]
        crds += [[conv(float(lines[2+i].split()[j])) for j in range(1,4)]]
    return (atoms, np.array(crds))


@engine.register_action
def read_xyzfile_crd(filename: "existing_file") -> "crd":
    conv = length_converter.get_converter('ang', 'au') 
    with open(filename, 'r') as infile:
        lines = infile.readlines()

    natoms = int(lines[0])
    crds = []
    for i in range(natoms):
        crds += [[conv(float(lines[2+i].split()[j])) for j in range(1,4)]]
    return np.array(crds)

@engine.register_action
def read_xyzfile_atomids(filename: "existing_file") -> "ilist":
    with open(filename, 'r') as infile:
        lines = infile.readlines()

    natoms = int(lines[0])
    atoms = []
    for i in range(natoms):
        atoms += [ATOMNAME_TO_ID[lines[2+i].split()[0]]]
    return atoms

