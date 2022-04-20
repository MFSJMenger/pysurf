from copy_execute import CopyExecute
from setup_spectrum import SetupSpectrum
from colt import Colt

class SubmitSpecCalc(Colt):
    """ Class to start single point calculation for the spectrum
        
        It uses the CopyExecute class and presets folder and subfolder
    """

    folder = SetupSpectrum.folder
    subfolder = SetupSpectrum.subfolder
    _user_input ="""
    copy = :: list
    exe = :: str
    """

    @classmethod
    def from_config(cls, config):
        config['folder'] = cls.folder
        config['subfolder'] = cls.subfolder
        return CopyExecute.from_config(config)

if __name__ == "__main__":
    SubmitSpecCalc.from_commandline()
