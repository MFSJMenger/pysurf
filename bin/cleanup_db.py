import numpy as np

from scipy.spatial.distance import cdist

from pysurf.colt import Colt
from pysurf.database import PySurfDB

class CleanupDB(Colt):
    _questions = """
    db_in = db.dat :: existing_file
    db_out = clean_db.dat :: file
    trust_radius_general = 0.75 :: float
    trust_radius_ci = 0.25 :: float
    #Energy difference in au, seperating CI trust radius and general trust radius
    energy_threshold = 0.02 :: float
    """

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        dbin = PySurfDB.load_database(config['db_in'], read_only=True)
        info = PySurfDB.info_database(config['db_in'])

        self.thresh = config['energy_threshold']
        self.trust_radius_general = config['trust_radius_general']
        self.trust_radius_ci = config['trust_radius_ci']

        if 'natoms' in info['dimensions']:
            model = False
        else:
            model = True
        dbout = PySurfDB.generate_database(config['db_out'], data=info['variables'], dimensions=info['dimensions'], model=model)

        self.crds = []
        for i, crd in enumerate(dbin['crd']):
            if len(dbout) != 0:
                diff = np.diff(dbin.get('energy', i))

                trust_general, trust_ci = self.is_within_radius(crd)
                if np.min(diff) < self.thresh:
                    if trust_ci is True:
                        continue
                else:
                    if trust_general is True:
                        continue

            self.crds += [crd]
            for prop in info['variables']:
                dbout.append(prop, dbin.get(prop, i))
            dbout.increase

    def is_within_radius(self, crd):
        crds = np.array(self.crds)
        shape = crds.shape
        if len(crds.shape) == 3:
            dist = cdist([np.array(crd).flatten()], crds.reshape((shape[0], shape[1]*shape[2])))
        else:
            dist = cdist([np.array(crd).flatten()], crds)
        trust_general = False
        trust_ci = False
        if np.min(dist) < self.trust_radius_general:
            trust_general = True
        if np.min(dist) < self.trust_radius_ci:
            trust_ci = True
        return trust_general, trust_ci

if __name__ == "__main__":
    CleanupDB.from_commandline()

