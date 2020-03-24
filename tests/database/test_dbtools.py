from pytest import fixture
import numpy as np

from pysurf.database.dbtools import DatabaseRepresentation


@fixture
def db_representation_basic():
    return """
      [dimensions]
      frames = unlimited
      natoms = 10
      """

@fixture
def db_representation_basic_changed():
    return """
      [dimensions]
      frames = unlimited
      natoms = 12
      """


def test_database_same_basic_representation(db_representation_basic): 
    rep1 = DatabaseRepresentation(db_representation_basic)
    rep2 = DatabaseRepresentation(db_representation_basic)
    assert(rep1 == rep2)


def test_database_different_basic_representation(db_representation_basic, db_representation_basic_changed): 
    rep1 = DatabaseRepresentation(db_representation_basic)
    rep2 = DatabaseRepresentation(db_representation_basic_changed)
    assert(rep1 != rep2)

@fixture
def db_representation_basic_changed():
    return {
      'dimensions': {
                'frames': 'unlimited', 
                'natoms': 12,
                }
    }
