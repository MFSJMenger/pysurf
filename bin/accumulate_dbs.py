import os
from shutil import copy2 as copy

from pysurf.colt import Colt
from pysurf.utils import SubfolderHandle
from pysurf.utils import exists_and_isfile
from pysurf.database import PySurfDB

from combine_dbs import CombineDBs

class AccumulateDBs(Colt):
    _questions = """
        #Foldername of the main folder, e.g. spectrum or prop
        folder = spectrum :: file

        #subfoldername, e.g. condition
        subfolder = condition :: str

        #db files that should be added to mother db
        dbfiles = db.dat :: str

        #mother db
        mother_db = db.dat :: str
        """


    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        setup = SubfolderHandle(config['folder'], config['subfolder'])

        counter = 0
        for file in setup.fileiter(config['dbfiles']):
            if counter == 0:
                copied = False
                if not exists_and_isfile(config['mother_db']):
                    copy(file, config['mother_db'])
                    copied = True
                info = PySurfDB.info_database(config['mother_db'])
                mother_db = PySurfDB.load_database(config['mother_db'], data=info['variables'], dimensions=info['dimensions'])
                counter += 1
                if copied is True:
                    continue

            CombineDBs(mother_db, file)
            print(f"Added file {file} to DB")
                    
                

if __name__=="__main__":
    AccumulateDBs.from_commandline()





