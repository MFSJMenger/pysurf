import configparser
import importlib
import os
# Numpy
import numpy as np
# Database related
from .dbinter.dbinter import DBInter
from .qminter.qminter import get_qminter
# utils
from ..utils.chemutils import atomic_masses
from ..utils.osutils import exists_and_isfile
from ..utils.context_utils import DoOnException
# fileparser
from ..fileparser import read_geom
# logger
from ..logger import get_logger, Logger


class SurfacePointProvider(object):
    """ The Surface Point Provider is the main interface providing the
        providing the information of a system at a specific point. It
        takes care where to take the information from according to the
        specified input file
    """

    def __init__(self, inputfile, logger=None):
        """ The inputfile for the SPP has to provide the necessary
            information, how to produce the data at a specific point
            in the crdinate space.
        """
        if not isinstance(logger, Logger):
            self.logger = get_logger('spp.log', 'SPP', [])
        else:
            self.logger = logger

        # get config
        self.config, self.path = self._parse_config(inputfile)
        # get current directory, which is the one where the actual calculation will be performed
        self.trajpath = os.getcwd()
        self.config['MAIN']['trajpath'] = self.trajpath
        #
        self.mode = self.config['MAIN'].get('mode', None)
        if self.mode is None:
            self.logger.error('Mode has to be provided in main section')
        #
        """ If a model is used, import the model according to the user
            input and provide an instance in the variable self.user_inst
        """
        if self.mode == 'model':
            self.logger.info('Using a model to generate the PES')
            self.interface = self._import_model_interface(self.config['MODEL']['module'],
                                                          self.config['MODEL']['class'])

        # If an ab initio calculation is used and if a database is used
        elif self.mode == 'ab initio':
            self.logger.info('Ab initio calculations are used to generate the PES')
            # make sure that AB INITIO section is in the inputfile
            # add path to AB INITIO section
            self.interface = self._import_abinitio_interface()
        else:
            self.logger.error("Mode has to be 'model' or 'ab initio'")

    def _parse_config(self, inputfile):
        """Parse the config file"""
        #
        config = configparser.ConfigParser()
        #
        if exists_and_isfile(inputfile):
            config.read(inputfile)
            inputfile_full = os.path.abspath(inputfile)
            path = os.path.dirname(inputfile_full)
            self.config['MAIN']['path'] = path
        else:
            self.logger.error('Inputfile '
                              + inputfile + ' for SurfacePointProvider '
                              + 'not found!')
        return config, path

    def get(self, request):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the crdinates and gives back the
            information at this specific position.
        """
        res = self.interface.get(request)
        # in the case of ab initio/DB add the masses
        if 'mass' in request and self.mode == 'ab initio':
            res['mass'] = self.get_masses()
        if 'atoms' in request and self.mode == 'ab initio':
            res['atoms'] = self.refgeo['atoms']
        return res

    def _get_refgeo(self, filename):
        atoms = []
        crds = []
        refgeo_path = os.path.join(self.path, filename)
        natoms, atoms, crds = read_geom(refgeo_path)
        return natoms, {'atoms': atoms, 'crd': np.array(crds)}

    def get_masses(self):
        masses = []
        for i in range(self.natoms):
            # masses are given in an array of shape (natoms, 3) like
            # the coordinates so that they can be easily used in the
            # surface hopping algorithm
            masses += [[atomic_masses[self.refgeo['atoms'][i]],
                        atomic_masses[self.refgeo['atoms'][i]],
                        atomic_masses[self.refgeo['atoms'][i]]]]
        return np.array(masses)

    def _error_on_exception(self, txt):
        """print error on exception and end code"""
        return DoOnException(self.logger.error, txt)

    def _import_model_interface(self, module_path, class_name):
        """Import Model class and initalize it"""
        path = os.abspath(module_path)
        self.logger.info('The model is given as: %s' % path)
        # load module
        with self._error_on_exception('The user module could not be loaded: %s' % path):
            spec = importlib.util.spec_from_file_location("", path)
            usr_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(usr_module)
        #
        with self._error_on_exception('The user class could not be found: %s' % class_name):
            usr_class = getattr(usr_module, class_name)
        #
        return usr_class()

    def _import_abinitio_interface(self):
        """Setup abinitio interface"""
        with self._error_on_exception('No AB INITIO section in the inputfile!'):
            self.config['AB INITIO']['path'] = self.path
            self.config['AB INITIO']['trajpath'] = self.trajpath
        # read reference geometry from inputfile
        refgeo_name = self.config['AB INITIO'].get('reference geometry', None)
        if refgeo_name is None:
            self.logger.error('no reference geometry file provided!')
        #
        self.natoms, self.refgeo = self._get_refgeo(refgeo_name)
        # get number of states
        if 'number of states' in self.config['AB INITIO']:
            with self._error_on_exception('Number of states is not an integer value!'):
                self.nstates = int(self.config['AB INITIO']['number of states'])
        else:
            self.logger.error('Number of states not specified in inputfile!')

        if 'database' in self.config['AB INITIO']:
            self.dbpath = os.path.join(self.path, self.config['AB INITIO']['database'])
            self.dbpath = os.path.realpath(self.dbpath)
            self.logger.info('Using database: ' + self.dbpath)
            self.db = True
            self.logger.add_handle("db")
            interface = DBInter(self.dbpath, self.config['AB INITIO'], self.logger["db"], self.refgeo)
        else:
            self.db = False
            self.logger.add_handle("QM")
            interface = get_qminter(self.config['AB INITIO'], self.logger["QM"], self.refgeo)
        return interface


if __name__ == "__main__":
    # spp = SurfacePointProvider('./test.inp')
    # print(spp.get(np.array([0.0, 0.0, 0.0])))
    spp = SurfacePointProvider('./test_abinit.inp')
    print(spp.get(spp.refgeo['ref geo']))
