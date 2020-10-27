import numpy as np

from scipy.spatial.distance import cdist
from scipy.spatial.distance import cdist, pdist, squareform

from colt import Colt
from pysurf.database import PySurfDB
from pysurf.spp import within_trust_radius
from pysurf.spp import internal


class CleanupDB(Colt):
    _questions = """
    db_in = db.dat :: file
    db_out = clean_db.dat :: file
    trust_radius_general = 0.75 :: float
    trust_radius_ci = 0.25 :: float
    #Energy difference in au, seperating CI trust radius and general trust radius
    energy_threshold = 0.02 :: float
    crd_mode = internal :: str :: [internal, cartesian]
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
        self.crd_mode = config['crd_mode']

        if 'natoms' in info['dimensions']:
            model = False
        else:
            model = True
        dbout = PySurfDB.generate_database(config['db_out'], data=info['variables'], dimensions=info['dimensions'], model=model)

        self.crds = None
        for i, crd in enumerate(dbin['crd']):
            if self.crd_mode == 'internal':
                crd = internal(np.copy(crd))
            else:
                crd = np.copy(crd)
            if i%1000 == 0:
                print(f"Processing point {i}")
            crd_shape = crd.shape
            crd.resize((1, crd.size))
            if self.crds is None:
                self.crds = crd.reshape((1, crd.size))
            else:
                diff = np.diff(dbin.get('energy', i))
                
                _, (trust_general, trust_ci) = within_trust_radius(crd, self.crds, radius=self.trust_radius_general, radius_ci=self.trust_radius_ci, metric='euclidean')
                if np.min(diff) < self.thresh:
                    if trust_ci is True:
                        continue
                else:
                    if trust_general is True:
                        continue
                self.crds = np.concatenate((self.crds, crd))

            crd.resize(crd_shape)
            for prop in info['variables']:
                dbout.append(prop, dbin.get(prop, i))
            dbout.increase

if __name__ == "__main__":
    CleanupDB.from_commandline()

