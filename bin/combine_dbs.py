
import numpy as np
from scipy.spatial import cKDTree

from pysurf.database import PySurfDB
from pysurf.colt import FromCommandline, Colt

""" This class has to be moved to spp.dbinter.dbinter.
    at the moment I cannot change dbinter, but it has to be done in futere
    so that the function can be just imported
"""
class NextNeighbor():
    def __init__(self, db):
        crds = []
        for crd in db['crd']:
            crds += [np.array(crd).flatten()]
        crds = np.array(crds)
        self.tree = cKDTree(crds)
    
    def get(self, crd):
        return self.tree.query(crd.flatten(), k=1)


class CombineDBs(Colt):
    _questions = """
    main_db = :: existing_file
    added_db = :: existing_file
    start_value = 0 :: int
    """

    @classmethod
    def from_config(cls, config):
        return cls(config['main_db'], config['added_db'], start=config['start_value'])

    
    def __init__(self, main_db, added_db, start=0):
    
        if not isinstance(main_db, PySurfDB):
            info = PySurfDB.info_database(main_db)
            if 'natoms' in info['dimensions']:
                model = False
            else:
                model = True
            main_db = PySurfDB.load_database(main_db, dimensions=info['dimensions'], data=info['variables'], model=model)
        if not isinstance(added_db, PySurfDB):
            added_db = PySurfDB.load_database(added_db, read_only=True)


        keys_raw = main_db.get_keys()
        keys = []
        for key in keys_raw:
            if main_db.get_dimension(key)[0].isunlimited() is True:
                keys += [key]

        #check whether dbs fit together
        check = True
        for key in keys:
            if key in added_db.get_keys():
                for dim1, dim2 in zip(main_db.get_dimension(key)[1:], added_db.get_dimension(key)[1:]):
                    if dim1.size != dim2.size:
                        check = False
            else:
                check = False
        if check is False:
            print('Error DBs do not fit together by dimensions')
            exit()
                


#        nn = NextNeighbor(db1)
        for i in range(start, len(added_db)):
#            min_dist = nn.get(crd)
#            if min_dist[0] > 0.25:
            for key in keys:
                main_db.append(key, added_db[key][i])
            main_db.increase
       

if __name__=='__main__':
    CombineDBs.from_commandline()
