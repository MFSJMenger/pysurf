import os
from shutil import copy2 as copy

from pysurf.colt import Colt
from pysurf.utils import SubfolderHandle

class CopyExecute(Colt):
    _questions = """
        #Foldername of the main folder, e.g. spectrum or prop
        folder = spectrum :: file

        #subfoldername, e.g. condition
        subfolder = condition :: str

        #list with files that should be copied
        copy = [submit.sh] :: list

        #executable that should be performed
        exe = sbatch submit.sh :: str
        """


    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        setup = SubfolderHandle(config['folder'], config['subfolder'])

        for subfolder in setup:
            print(subfolder)
            for item in config['copy']:
                copy(item, subfolder)
        
            os.chdir(subfolder)
            os.system(config['exe'])
            os.chdir(setup.main_folder)

if __name__=="__main__":
    CopyExecute.from_commandline()





