from pytest import fixture
import os
import numpy as np

from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable, DatabaseRepresentation


@fixture
def filepath():
    return 'database.nc'

@fixture
def default_settings():
    return {
      'dimensions': {
                'frames': 'unlimited', 
                'natoms': 10,
                'three': 3,
                'one': 1,
       },
       'variables': {
                'dipole': DBVariable(np.double, ('frames', 'three'))
       }
    }


@fixture
def empty_database(filepath, default_settings):

    if os.path.exists(filepath):
        if os.path.isfile(filepath):
            os.remove(filepath)
        else:
            raise Exception(f"file '{filepath}' exists is not a file")

    return Database(filepath, default_settings)


@fixture
def empty_database_closed(empty_database):
    empty_database.close()
    return empty_database


def test_check(empty_database, default_settings):
    rep = DatabaseRepresentation.from_db(empty_database._db)
    rep2 = DatabaseRepresentation(default_settings)
    assert(rep == rep2)
    
def test_append(empty_database_closed, filepath, default_settings):
    nc = Database(filepath, default_settings)
    nc.close()
    nc = Database(filepath, default_settings)
    nc.append('dipole', np.array([1, 2, 3]))
    assert(nc.get_dimension_size('frames') == 1)


    
