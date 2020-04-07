import os
# Numpy
import numpy as np
# Database related
from ..colt import Colt, PluginBase
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

"""


class NoFurtherQuestions(Colt):
    _questions = ""


class AbinitioFactory(PluginBase):
    """Factory for any QM code"""

    _is_plugin_factory = True
    _plugins_storage = 'software'

    _questions = """
        software =
    """

    reader = None

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.add_branching("software", {name: software.questions for name, software in cls.software.items()})

    @classmethod
    def instance_from_config(cls, config, atomids, nstates):
        return cls.software[config['software'].value].from_config(config['software'].subquestion_answers, atomids, nstates)


class ModelFactory(PluginBase):
    """Model Factory"""

    _is_plugin_factory = True
    _plugins_storage = '_models'

    _questions = """
        model =
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("model", {name: model.questions for name, model in cls._models.items()})

    @classmethod
    def instance_from_config(cls, config):
        model = config['model'].value
        return cls._models[model].from_config(config['model'])


class SurfacePointProvider(Colt):
    """ The Surface Point Provider is the main interface providing the
        providing the information of a system at a specific point. It
        takes care where to take the information from according to the
        specified input file
    """

    _questions = """
        logging = debug :: str ::
        mode = ab-initio
        use_db = no 
        """
    _modes = {'ab-initio': AbinitioFactory,
              'model': ModelFactory,
    }

    _database = {
            'yes': DataBaseInterpolation,
            'no': NoFurtherQuestions
            }

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("mode", {name: mode.questions for name, mode in cls._modes.items()})
        questions.generate_cases("use_db", {name: mode.questions for name, mode in cls._database.items()})

    def __init__(self, inputfile, properties, natoms, nstates, atomids, logger=None):
        """ The inputfile for the SPP has to provide the necessary
            information, how to produce the data at a specific point
            in the coordinate space.

            Args:

                inputfile, str:
                    Name of the input file for the SPP

                properties, list/tuple
                    all possible requested properties, used for sanity checks

                natoms, int:
                     Number of atoms in the system

                nstates, int:
                    Number of states requested

                logger, optional
                    Logging module used
        """

        if not isinstance(logger, Logger):
            self.logger = get_logger('spp.log', 'SPP', [])
        else:
            self.logger = logger

        # get config
        config = self._parse_config(inputfile)
        #
        self._request, self._interface = self._select_interface(config, properties, natoms, nstates, atomids)

    def _select_interface(self, config, properties, natoms, nstates, atomids):
        """Select the correct interface based on the mode"""
        if config['mode'] == 'model':
            self.logger.info('Using model to generate the PES')
            interface = ModelFactory.instance_from_config(config['mode'])
        # If an ab initio calculation is used and if a database is used
        elif config['mode'] == 'ab-initio':
            self.logger.info('Ab initio calculations are used to generate the PES')
            # make sure that AB INITIO section is in the inputfile
            # add path to AB INITIO section
            interface = AbinitioFactory.instance_from_config(config['mode'], atomids, nstates)
        else:
            # is atomatically caught through colt!
            self.logger.error("Mode has to be 'model' or 'ab-initio'")
        # check default
        self._check_properties(properties, interface)
        # use databse
        if config['use_db'] == 'yes':
            self.logger.info("Setting up database...")
            interface = DataBaseInterpolation(interface, config['use_db'], natoms, nstates, properties)
            request = RequestGenerator(nstates, properties + config['use_db']['properties'])
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

    def _parse_config(self, inputfile):
        """Parse the config file"""
        questions = self.generate_questions("spp", config=None)
        return questions.check_only(inputfile)

    def request(self, crd, properties, states=None):
        request = self._request.request(crd, properties, states)
        return self.get(request)

    def get(self, request):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the crdinates and gives back the
            information at this specific position.

            It does not perform any sanity checks anylonger, so insure that all
            possible requested properties are in used!
        """
        return self._interface.get(request)
