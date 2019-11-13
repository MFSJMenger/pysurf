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
