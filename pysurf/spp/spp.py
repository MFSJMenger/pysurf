import os
# Numpy
import numpy as np
# Database related
from colt import Colt, Plugin
from colt.obj import NoFurtherQuestions
# logger
from ..logger import get_logger, Logger
# Interpolation
from .dbinter import DataBaseInterpolation
#
from .request import RequestGenerator

"""
TODO:
    - What happens to user properties that are asked
      if no database is set?
    - Models should check that enough states are computed!
"""


class AbinitioFactory(Plugin):
    """Factory for any QM code"""

    _is_plugin_factory = True
    _plugins_storage = 'software'

    _user_input = """
        software =
    """

    reader = None

    @classmethod
    def _extend_user_input(cls, questions):
        questions.add_branching("software", {name: software.colt_user_input
                                             for name, software in cls.software.items()})


class ModelFactory(Plugin):
    """Model Factory"""

    _is_plugin_factory = True
    _plugins_storage = '_models'

    _user_input = """
        model =
    """

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases("model", {name: model.colt_user_input
                                           for name, model in cls._models.items()})


class SurfacePointProvider(Colt):
    """ The Surface Point Provider is the main interface providing
        to perform calculations of a given system at a specific coordinate.
    """
    # Main Questions for SPP
    _user_input = """
        logging = debug :: str ::
        mode = ab-initio
        use_db = no
        """
    # different modes
    _modes = {'ab-initio': AbinitioFactory,
              'model': ModelFactory,
    }
    # use_db
    _database = {
            'yes': DataBaseInterpolation,
            'no': NoFurtherQuestions
            }

    @classmethod
    def _extend_user_input(cls, questions):
        """ Extend the questions of the surface point provider using
            the questions of the `Abinitio` and `Model` """
        questions.generate_cases("mode", {name: mode.colt_user_input for name, mode in cls._modes.items()})
        questions.generate_cases("use_db", {name: mode.colt_user_input for name, mode in cls._database.items()})

    @classmethod
    def from_config(cls, config, properties, nstates, natoms, *,
                    nghost_states=0, atomids=None, logger=None):
        """ Initialize the SurfacePointProvider from the information
            read from configfile using `Colt`

            Parameters
            ----------

            config: ColtConfig
                config object created through colt by reading a configfile/asking
                commandline questions, etc.

            properties: list/tuple of strings
                all possible requested properties, used for sanity checks

            natoms: int
                Number of atoms in the system

            nstates: int
                Number of states requested

            nghost_states: int
                Number of ghost states requested

            atomids: list/tuple, optional
                Atomids of the molecule, needed in case an abinitio method is selected
                ignored in case of the model

            logger: Logger, optional
                objected used for loggin, if None, a new one will be created

            Returns
            -------

            SurfacePointProvider
                SurfacePointProvider instance

        """

        return cls(config['mode'], config['use_db'],
                   properties, nstates, natoms,
                   nghost_states=nghost_states, atomids=atomids, logger=logger,
                   logging_level=config['logging'])

    def __init__(self, mode_config, use_db, properties, nstates, natoms, *,
                 nghost_states=0, atomids=None, logger=None, logging_level='debug'):
        """ The inputfile for the SPP has to provide the necessary
            information, how to produce the data at a specific point
            in the coordinate space.

            Parameters
            ----------

                mode_config:
                    `config` stating the interface etc. and corresponding info

                use_db:
                    `config` stating if database should be setup or not

                properties: list/tuple of strings
                    all possible requested properties, used for sanity checks

                natoms: int
                     Number of atoms in the system

                nstates: int
                    Number of states requested

                nghost_states: int
                    Number of ghost states requested

                atomids: list/tuple, optional
                    Atomids of the molecule, needed in case an abinitio method is selected
                    ignored in case of the model

                logging_level: str, optional
                    Specifies the logging level of the `SurfacePointProvider`

                logger: Logger, optional
                    objected used for loggin, if None, a new one will be created
        """
        if logger is None:
            self.logger = get_logger('spp.log', 'SPP', [])
        else:
            self.logger = logger
        #
        self._request, self._interface = self._select_interface(mode_config, use_db, 
                                                                properties, natoms, 
                                                                nstates, nghost_states, 
                                                                atomids)

    @property
    def interpolator(self):
        interpolator = getattr(self._interface, 'interpolator', None)
        if interpolator is None:
            raise Exception("Interpolator not available")
        return interpolator

    def request(self, crd, properties, states=None, same_crd=False):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the crdinates and gives back the
            information at this specific position.

            It does not perform any sanity checks anylonger, so insure that all
            possible requested properties are in used!
        """
        datacontainer = self._request.request(crd, properties, states, same_crd=same_crd)
        return self._interface.get(datacontainer)

    def _select_interface(self, mode_config, use_db, properties, natoms, 
                          nstates, nghost_states, atomids):
        """Select the correct interface based on the mode"""
        #
        if mode_config == 'model':
            self.logger.info('Using model to generate the PES')
            interface = ModelFactory.plugin_from_config(mode_config['model'])
        # If an ab initio calculation is used and if a database is used
        elif mode_config == 'ab-initio':
            self.logger.info('Ab initio calculations are used to generate the PES')
            # make sure that AB INITIO section is in the inputfile
            # add path to AB INITIO section
            interface = AbinitioFactory.plugin_from_config(mode_config['software'], atomids, nstates, nghost_states)
        else:
            # is atomatically caught through colt!
            self.logger.error("Mode has to be 'model' or 'ab-initio'")
        # check default
        self._check_properties(properties, interface)
        # use databse
        if use_db == 'yes':
            self.logger.info("Setting up database...")
            interface = DataBaseInterpolation.from_config(use_db, interface, natoms, nstates, 
                                              properties, model=(mode_config=='model'))
            if use_db['properties'] is not None:
                properties += use_db['properties']
            request = RequestGenerator(nstates, properties, use_db=True)
            self.logger.info("Database ready to use")
        else:
            request = RequestGenerator(nstates)
        #
        return request, interface

    def _check_properties(self, properties, interface):
        """Sanity check for properties"""
        if any(prop not in interface.implemented for prop in properties):
            self.logger.error(f"""The Interface does not provide all properties:
Needed: {properties}
Implemented: {interface.implemented}
            """)
        self.logger.info(f"Interface {interface} provides all necessary properties")


