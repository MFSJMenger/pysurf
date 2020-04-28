from copy_execute import CopyExecute
from setup_propagation import SetupPropagation
from pysurf.colt import Colt

class SubmitPropCalc(Colt):
    """ Class to start single point calculation for the spectrum
        
        It uses the CopyExecute class and presets folder and subfolder
    """

    folder = SetupPropagation.folder
    subfolder = SetupPropagation.subfolder
    _questions="""
    copy = :: list
    exe = :: str
    """

    @classmethod
    def from_config(cls, config):
        config['folder'] = cls.folder
        config['subfolder'] = cls.subfolder
        return CopyExecute.from_config(config)

if __name__ == "__main__":
    SubmitPropCalc.from_commandline()
