from shutil import copy2 as copy
from subprocess import run, CalledProcessError

from colt import Colt
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
            for item in config['copy']:
                copy(item, subfolder)
        
            try:
                run(config['exe'], cwd=subfolder, check=True, shell=True)
            except KeyboardInterrupt or CalledProcessError:
                break


if __name__=="__main__":
    CopyExecute.from_commandline()
