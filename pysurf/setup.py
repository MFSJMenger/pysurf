from abc import abstractmethod
#
from .utils import SubfolderHandle
from .logger import get_logger
from .colt import Colt


class SetupBase(Colt, SubfolderHandle):

    folder = None
    subfolder = None

    def __init__(self, logger=None, digit=8):
        if self.folder is None:
            raise Exception("self.folder needs to be set!")
        if self.subfolder is None:
            raise Exception("self.subfolder needs to be set!")
        #
        SubfolderHandle.__init__(self, self.folder, self.subfolder)
        #
        if logger is None:
            self.logger = get_logger(None, "")
        else:
            self.logger = logger

    def setup_folders(self, lst, *args, **kwargs):
        self.logger.info(f'Start setting up folders\n')
        for i, folder in enumerate(self.folderiter(lst)):
            self.logger.info(f'Create folder {folder}\n')
            self.setup_folder(i, folder, *args, **kwargs)

    @abstractmethod                
    def setup_folder(self, number, foldername, *args, **kwargs):
        """fill the folder with data"""
