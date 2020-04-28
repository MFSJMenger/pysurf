import os
import re
from .osutils import exists_and_isfile

def create_folder(folder): 
    if os.path.exists(folder):
        if not os.path.isdir(folder):
            raise Exception(f"{folder} needs to be a folder")
    else:
        os.mkdir(folder)


class SubfolderHandle:
    """Creates a useful handle to an folder tree of following structure
       
       folder/
             subfolder_0000001/
                               ouput.txt
                               input.txt
                               value.txt
             subfolder_0000002/
             subfolder_0000003/
             ...
             subfolder_0001000/
             subfolder_0001001/
             subfolder_0001002/
             ...
        
    """

    def __init__(self, folder, subfolder, digits=8): 
        """ """
        self.parent = os.getcwd()
        self.main_folder = self._create_main_folder(folder)
        self.template, self.reg = self._setup(subfolder, digits)
        self._sanity_check()
        self._folders = self._get_folders()

    def fileiter(self, filename): 
        """Returns an iterator over all files with name
           filename in folder/subfolder_*/ """
        for folder in self:
            name = os.path.join(folder, filename)
            if os.path.isfile(name):
                yield name
                
    def setupiter(self, lst): 
        """Returns an iterator over all files with name
           filename in folder/subfolder_*/ """
        for idx in lst:
            name = self._folder_path(idx)
            if not os.path.exists(name):
                create_folder(name)
                yield idx, name

    def folderiter(self, lst):
        """iterates of the folders idx and names of folders 
           with a certain idx which do not exists yet and creates them
         """
        for idx in lst:
            name = self._folder_path(idx)
            create_folder(name)
            yield name
        self._update_folders()

    def get_file(self, filename, idx):
        filepath = os.path.join(self.main_folder, self._folder_name(idx), filename)
        if exists_and_isfile(filepath):
            return filepath
        else:
            return None

    def generate_folders(self, lst):
        """Generates folders"""
        for idx in lst:
            create_folder(self._folder_path(idx))
        self._update_folders()

    def __iter__(self):
        for folder in self._folders:
            yield folder

    def __len__(self):
        return len(self._folders)

    def _folder_path(self, idx):
        """Return absolute path to an given folder"""
        return os.path.join(self.main_folder, self._folder_name(idx))

    def _folder_name(self, idx):
        return self.template % idx

    def _update_folders(self):
        self._folders = self._get_folders()

    def _setup(self, subfolder, digits):
        template = f"{subfolder}_%0{digits}d"
        # we except any digits in case one changes those
        regex = f"{subfolder}_" + r"\d+"
        reg = re.compile(regex)
        return template, reg

    def _create_main_folder(self, folder):
        folder = os.path.abspath(folder) 
        create_folder(folder)
        return folder

    def _get_folders(self): 
        return tuple(sorted(os.path.join(self.main_folder, filename)
                            for filename in os.listdir(self.main_folder)
                            if self.reg.match(filename) is not None))

    def _sanity_check(self):
        name = self._folder_name(1)
        if self.reg.match(name) is None:
            raise Exception("template and regex are not in sync")
