class Molecule:
    """Store info of a molecule"""
    def __init__(self, atomids, crd, masses=None, name=None):
        self.natoms = len(atomids)
        self.atomids = atomids
        self.crd = crd
        self.name = name
        self.masses = masses

    def __len__(self):
        return self.natoms

    def __iter__(self):
        return zip(self.atomids, self.crd)

    def write_xyz(self, filename):
        natoms = len(self.atomids)
        out = f"""{natoms}
        commentar\n"""
        out += "\n".join(self.format(atomid, crd) for atomid, crd in self)
        with open(filename, 'w') as f:
            f.write(out)

    def format(self, atomid, crd):                    
        return "%s   %12.8f %12.8f %12.8f" % (atomid, crd[0], crd[1], crd[2])
