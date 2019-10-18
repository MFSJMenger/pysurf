import configparser
import os
import numpy as np
import importlib
import logging

from .qminter.qminter import get_qminter
from ..utils.chemutils import atomic_masses
from ..database.database import Database
from ..database.dbtools import DBVariable
from .dbinter.dbinter import DBInter
from ..utils.constants import angstrom2bohr

class SurfacePointProvider(object):
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

        #close all open handlers of spp loggers
        for hand in self.logger.handlers:
           hand.stream.close()
           self.logger.removeHandler(hand)

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

        # get current directory, which is the one where the abinit
        # calculation will be performed
        self.trajpath = os.getcwd()
        self.config['MAIN']['trajpath'] = self.trajpath


        if 'logging' in self.config['MAIN'].keys():
            loglevel = self.config['MAIN']['logging']
            if loglevel == 'debug':
                #handler.setLevel(logging.DEBUG)
                self.logger.setLevel(logging.DEBUG)
                for hand in self.logger.handlers:
                    hand.setLevel(logging.DEBUG)
                self.logger.info('logging level set to debug')
            else:
                self.logger.setLevel(logging.INFO)
                self.logger.info('logging level set to info' +
                                 ' as specification was not known')

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
            self.interface = usr_class()

        #If an ab initio calculation is used
        #    and if a database is used
        elif self.mode == 'ab initio':
            self.logger.info('Ab initio calculations are used to '
                             + 'generate the PES')

            # make sure that AB INITIO section is in the inputfile
            # add path to AB INITIO section
            try:
                self.config['AB INITIO']['path'] = self.path
                self.config['AB INITIO']['trajpath'] = self.trajpath

            except KeyError:
                self.logger.error('No AB INITIO section in '
                                  + 'the inputfile!')
                exit()

            # read reference geometry from inputfile
            self.refgeo = self.get_refgeo()
            self.natoms = len(self.refgeo['atoms'])

            if 'number of states' in self.config['AB INITIO'].keys():
                try:
                    self.nstates = int(self.config['AB INITIO']
                                       ['number of states'])
                except ValueError:
                    self.logger.error('Number of states is not an'
                                      + 'integer value!')
                    exit()
            else:
                self.logger.error('Number of states not specified'
                                   + 'in inputfile!')
                exit()
            if 'database' in self.config['AB INITIO'].keys():
                self.dbpath = os.path.join(self.path, self.config['AB INITIO']['database'])
                self.dbpath = os.path.realpath(self.dbpath)
                self.logger.info('Using database: '
                                 + self.dbpath)
                self.interface = DBInter(self.dbpath, self.config['AB INITIO'], self.logger, self.refgeo)
                self.db = True
            else:
                self.db = False

                self.interface = get_qminter(self.config['AB INITIO'],
                                      self.logger, self.refgeo)

    def get(self, request):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the coordinates and gives back the
            information at this specific position.
        """
        res = self.interface.get(request)
        # in the case of ab initio/DB add the masses if not done by
        # the interface
        if 'mass' in request.keys() and self.mode == 'ab initio':
            res['mass'] = self.get_masses()

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
                    coords += [[float(c)*angstrom2bohr for c in split_line[1:]]]
        for i in range(len(atoms)):
            atom = atoms[i]
            atom = atom[0].upper() + atom[1:].lower()
            atoms[i] = atom

        return {'atoms': atoms, 'coord': np.array(coords)}

    def get_masses(self):
        masses = []
        for i in range(len(self.refgeo['atoms'])):
            # masses are given in an array of shape (natoms, 3) like
            # the coordinates so that they can be easily used in the
            # surface hopping algorithm
            masses += [[atomic_masses[self.refgeo['atoms'][i]],
                        atomic_masses[self.refgeo['atoms'][i]],
                        atomic_masses[self.refgeo['atoms'][i]]]]
        return np.array(masses)


if __name__ == "__main__":
    # spp = SurfacePointProvider('./test.inp')
    # print(spp.get(np.array([0.0, 0.0, 0.0])))
    spp = SurfacePointProvider('./test_abinit.inp')
    print(spp.get(spp.refgeo['ref geo']))
