#! /data/ehrmaier/anaconda3/bin/python3
import os
import numpy as np

from pysurf.database.database import Database
from pysurf.colt import FromCommandline

@FromCommandline("""
infile = :: file
key = :: str
""")
def open_db_command(infile, key):
    open_db(infile, key)

def open_db(infile, key):    
    if os.path.isfile(infile):
        pass
    else:
        print('Error: Cannot find database file!')
        exit()
    
    db = Database.load_db(infile)
    
    
    #for crd in db['crd']:
    #    print(crd)
    
    #for veloc in db['veloc']:
    #    print(veloc)
    keys = db.get_keys()
    print('keys in the database: ')
    string = ''
    for k in keys:
        string += ' ' + k
    print(string)
    print('entries in the db: ', len(db[k]))

    print('printing all entries of '+key+':')
    for value in db[key]:
        print(value)

if __name__=="__main__":
    open_db_command()
