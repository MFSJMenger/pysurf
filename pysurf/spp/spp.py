import configparser
import os
import numpy as np
import importlib
import logging

from .interface.interface import Interface
from ..utils.chemutils import atomic_masses

class SurfacePointProvider():
    """ The Surface Point Provider is the main interface providing the
        providing the information of a system at a specific point. It
        takes care where to take the information from according to the
        specified input file
    """
    def __init__(self, inputfile):
        """ The inputfile for the SPP has to provide the necessary
            information, how to produce the data at a specific point
            in the coordinate space.
        """

        self.logger = logging.getLogger('spp')
        self.logger.setLevel(logging.INFO)
        # create a file handler
        handler = logging.FileHandler(filename='spp.log', mode='w')
        handler.setLevel(logging.INFO)
        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(name)s - '
                                      + '%(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        # add the file handler to the logger
        self.logger.addHandler(handler)
        self.logger.info('Surface Point Provider')

        self.config = configparser.ConfigParser()
        if os.path.isfile(inputfile):
            self.config.read(inputfile)
            inputfile_full = os.path.abspath(inputfile)
            self.path = os.path.dirname(inputfile_full)
            self.config['MAIN']['path'] = self.path
        else:
            self.logger.error('Inputfile '
                              + inputfile + ' for SurfacePointProvider '
                              + 'not found!')
            exit()
        if 'mode' in self.config['MAIN'].keys():
            self.mode = self.config['MAIN']['mode']
        else:
            self.logger.error('Mode has to be provided in main section')
            exit()

        """ If a model is used, import the model according to the user
            input and provide an instance in the variable self.user_inst
        """
        if self.mode == 'model':
            self.logger.info('Using a model to generate the PES')

            try:
                # get absolute path for model
                path = os.path.join(self.path, self.config['MODEL']['module'])
                path = os.path.realpath(path)
                self.logger.info('The model is given as: ' + path)
                spec = importlib.util.spec_from_file_location("", path)
                usr_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(usr_module)
            except (FileNotFoundError, ModuleNotFoundError, ImportError):
                self.logger.error('The user module could not be '
                                  + 'loaded: ' + path)
                exit()
            try:
                usr_class = getattr(usr_module, self.config['MODEL']['class'])
            except AttributeError:
                self.logger.error('The user class could not be found: '
                                  + self.config['MODEL']['class'])
                exit()
            self.usr_inst = usr_class()
        elif self.mode == 'ab initio':
            self.logger.info('Ab initio calculations are used to '
                             + 'generate the PES')
            # read reference geometry from inputfile
            self.refgeo = self.get_refgeo()
            if 'database' in self.config['MAIN'].keys():
                self.db = self.config['MAIN']['database']
                self.logger.info('Using database: '
                                 + self.db)
            else:
                self.db = False

    def get(self, coord):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the coordinates and gives back the
            information at this specific position.
        """
        if self.mode == 'model':
            return self.usr_inst.get(coord)
        elif self.mode == 'ab initio':
            if self.db is False:
                return self.qm_get(coord)
            else:
                self.logger.error('DB not yet implemented!')
                exit()

    def qm_get(self, coord):
        try:
            self.config['AB INITIO']['path'] = self.path
            interface = Interface(self.config['AB INITIO'],
                                  self.logger, self.refgeo)
        except KeyError:
            self.logger.error('No AB INITIO section in '
                              + 'the inputfile!')
            exit()
        res = interface.get(coord)
        masses = self.get_masses()
        res['mass'] = masses
        return res

    def get_refgeo(self):
        if 'reference geometry' not in self.config['AB INITIO'].keys():
            self.logger.error('no reference geometry file provided!')
            exit()
        atoms = []
        coords = []
        refgeo_path = os.path.join(self.path,
                                   self.config['AB INITIO']
                                   ['reference geometry'])
        with open(refgeo_path) as infile:
            infile.readline()
            infile.readline()
            for line in infile:
                split_line = line.split()
                if len(split_line) == 4:
                    atoms += [split_line[0]]
                    coords += [[float(c) for c in split_line[1:]]]
        for i in range(len(atoms)):
            atom = atoms[i]
            atom = atom[0].upper() + atom[1:].lower()
            atoms[i] = atom

        return {'atoms': atoms, 'coord': np.array(coords)}

    def get_masses(self):
        masses = []
        for i in range(len(self.refgeo['atoms'])):
            masses += [atomic_masses[self.refgeo['atoms'][i]]]
        return np.array(masses)


if __name__ == "__main__":
    # spp = SurfacePointProvider('./test.inp')
    # print(spp.get(np.array([0.0, 0.0, 0.0])))
    spp = SurfacePointProvider('./test_abinit.inp')
    print(spp.get(spp.refgeo['ref geo']))
