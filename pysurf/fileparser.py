from .atominfo import atomid_to_name
#
from .utils.osutils import exists_and_isfile
from .utils.constants import angstrom2bohr


def read_geom(filename):

    if not exists_and_isfile(filename):
        raise Exception("File '%s' does not exisits" % filename)

    ftype = filename.rpartition('.')[2]

    if ftype == 'xyz':
        return _read_xyz(filename)

    raise Exception("Filetype '%s' unkown!" % ftype)


def _read_xyz(filename):
    atoms = []
    coords = []
    with open(filename, 'r') as f:
            natoms = int(f.readline())
            f.readline()
            for _ in range(natoms):
                line = f.readline()
                split_line = line.split()
                if len(split_line) == 4:
                    atoms.append(split_line[0])
                    coords.append([float(c)*angstrom2bohr for c in split_line[1:5]])

    atoms_are_numbers = False

    try:
        a = int(atoms[0])
        atoms_are_numbers = True
    except:
        pass

    if atoms_are_numbers is True:
        for idx, atom in enumerate(atoms):
            atoms[idx] = atomid_to_name[int(atom)
    else:
        for idx, atom in enumerate(atoms):
            atom = atom[0].upper() + atom[1:].lower()
            atoms[idx] = atom

    return natoms, atoms, coords
