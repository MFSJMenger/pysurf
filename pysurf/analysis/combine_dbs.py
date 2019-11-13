import click
import numpy as np
from scipy.spatial import cKDTree

from pysurf.database.database import Database
#from pysurf.spp.dbinter.dbinter import get_min_dist

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

@click.command()
@click.argument('database1')
@click.argument('database2')
def combine_dbs_command(database1, database2):
    combine_dbs(database1, database2)

def combine_dbs(database1, database2):
    print('database 2 is added to database 1')
    print('database1 ', database1)
    print('database2 ', database2)

    db1 = Database.load_db(database1)
    db2 = Database.load_db(database2)
    thr = 0.25

    keys = db2.get_keys()
    nn = NextNeighbor(db1)
    for i,crd in enumerate(db2['crd']):
        min_dist = nn.get(crd)
        if min_dist[0] > 0.25:
            for key in keys:
                db1.append(key,db2[key][i])
            db1.increase
       

if __name__=='__main__':
    combine_dbs_command()
