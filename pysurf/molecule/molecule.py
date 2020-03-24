class Molecule(object):
    """Store info of a molecule"""
    def __init__(self, atomids, crd, masses, name=None):
        self.atomids = atomids
        self.crd = crd
        self.masses = masses
        self.name = name

    @property
    def natoms(self):
        return len(self.atomids)

    def write_xyz(self, filename):
        natoms = len(self.atomids)
        out = f"""{natoms}
        commentar\n"""
        out += "\n".join(self.format(idx) for idx in range(natoms))
        with open(filename, 'w') as f:
            f.write(out)

    def format(self, idx):                    
        return "%s   %12.8f %12.8f %12.8f" % (self.atomids[idx], self.crd[idx][0], self.crd[idx][1], self.crd[idx][2])
    
