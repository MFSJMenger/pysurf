import numpy as np

from scipy.spatial.distance import cdist

from pysurf.colt import Colt
from pysurf.database import PySurfDB

class CleanupDB(Colt):
    _questions = """
    db_in = db.dat :: existing_file
    db_out = clean_db.dat :: file
    threshold = 0.25 :: float
    """

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        dbin = PySurfDB.load_database(config['db_in'], read_only=True)
        info = PySurfDB.info_database(config['db_in'])

        self.thresh = config['threshold']

        if 'natoms' in info['dimensions']:
            model = False
        else:
            model = True
        dbout = PySurfDB.generate_database(config['db_out'], data=info['variables'], dimensions=info['dimensions'], model=model)

        self.crds = []
        for i, crd in enumerate(dbin['crd']):
            if i != 0:
                if self.is_within_radius(crd) is True:
                    continue
            self.crds += [crd]
            for prop in info['variables']:
                dbout.append(prop, dbin.get(prop, i))
            dbout.increase

    def is_within_radius(self, crd):
        dist = cdist([crd], self.crds)
        if np.min(dist) < self.thresh:
            return True
        else:
            return False

if __name__ == "__main__":
    CleanupDB.from_commandline()

