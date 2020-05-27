import numpy as np

from scipy.spatial.distance import cdist

from pysurf.colt import Colt
from pysurf.database import PySurfDB
from scipy.spatial.distance import cdist, pdist, squareform


def inverse(crd):
    return pdist(crd)
    return np.array([1.0/ele for ele in pdist(crd)])


def inverse_coordinates(crds):
    return np.array([inverse(crd) for crd in crds])


class CleanupDB(Colt):
    _questions = """
    db_in = db.dat :: existing_file
    db_out = clean_db.dat :: file
    trust_radius_general = 0.75 :: float
    trust_radius_ci = 0.25 :: float
    #Energy difference in au, seperating CI trust radius and general trust radius
    energy_threshold = 0.02 :: float
    inverse = false :: bool
    """

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        dbin = PySurfDB.load_database(config['db_in'], read_only=True)
        info = PySurfDB.info_database(config['db_in'])

        self.inverse = config['inverse']
        self.thresh = config['energy_threshold']
        self.trust_radius_general = config['trust_radius_general']
        self.trust_radius_ci = config['trust_radius_ci']

        if 'natoms' in info['dimensions']:
            model = False
        else:
            model = True
        dbout = PySurfDB.generate_database(config['db_out'], data=info['variables'], dimensions=info['dimensions'], model=model)

        self.crds = None
        for i, crd in enumerate(dbin['crd']):
            if i%1000 == 0:
                print(f"Processing point {i}")
            crd = np.copy(crd)
            crd_shape = crd.shape
            crd.resize((1, crd.size))
            if self.crds is None:
                self.crds = crd.reshape((1, crd.size))
            else:
                diff = np.diff(dbin.get('energy', i))
                
                trust_general, trust_ci = self.is_within_radius(crd)
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

    def is_within_radius(self, crd):
        if self.inverse is True:
            crd = inverse(crd)
            crds = inverse_coordinates(np.array(self.crds))
        else:
            crds = self.crds
        #
        shape = crds.shape
        shape_crd = crd.shape
        if len(crds.shape) == 3:
            dist = cdist(crd.resize((1, crd.size)), crds.resize((shape[0], shape[1]*shape[2])), metric=dim_norm)
        else:
            crd.resize((1, crd.size))
            dist = cdist(crd, crds, metric=dim_norm)
        crd.resize(shape_crd)
        trust_general = False
        trust_ci = False
        if np.min(dist) < self.trust_radius_general:
            trust_general = True
        if np.min(dist) < self.trust_radius_ci:
            trust_ci = True
        return trust_general, trust_ci


def dim_norm(crd1, crd2):
    return np.max(np.abs(crd1-crd2))

if __name__ == "__main__":
    CleanupDB.from_commandline()

