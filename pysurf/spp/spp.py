import os
# Numpy
import numpy as np
# Database related
#
from ..colt import Colt, PluginBase
from ..molecule.atominfo import masses_au as atomic_masses
# fileparser
from ..fileparser import read_geom
# logger
from ..logger import get_logger, Logger


class AbinitioFactory(PluginBase):
    """Factory for any QM code"""

    _is_plugin_factory = True
    _plugins_storage = 'software'

    _questions = """
        software = 
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.add_branching("software", {name: software.questions for name, software in cls.software.items()})

    @classmethod
    def class_from_config(cls, config):
        return cls.software[config['software'].value](config['software'].subquestion_answers)


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
    def class_from_config(cls, config):
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

        [use_db(yes)]
        properties = NO :: str

        [use_db(no)]
        """
    _modes = {'ab-initio': AbinitioFactory,
              'model': ModelFactory, 
    }

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases("mode", {name: mode.questions for name, mode in cls._modes.items()})

    def __init__(self, inputfile, logger=None):
        """ The inputfile for the SPP has to provide the necessary
            information, how to produce the data at a specific point
            in the coordinate space.
        """

        if not isinstance(logger, Logger):
            self.logger = get_logger('spp.log', 'SPP', [])
        else:
            self.logger = logger

        # get config
        config = self._parse_config(inputfile)
        #
        self._interface = self._select_interface(config)

    def _select_interface(self, config):

        if config['mode'] == 'model':
            self.logger.info('Using model to generate the PES')
            interface = ModelFactory.class_from_config(config['mode'])
        # If an ab initio calculation is used and if a database is used
        elif config['mode'] == 'ab-initio':
            self.logger.info('Ab initio calculations are used to generate the PES')
            # make sure that AB INITIO section is in the inputfile
            # add path to AB INITIO section
            interface = AbinitioFactory.class_from_config(config['mode'])
        else:
            self.logger.error("Mode has to be 'model' or 'ab-initio'")
        return interface

    def _parse_config(self, inputfile):
        """Parse the config file"""
        with self.logger.exit_on_exception(f"Inputfile '{inputfile} "
                                           "for SurfacePointProvider not found!"):
            questions = self.generate_questions("spp", inputfile)
        #
        return questions.ask()

    def get(self, request):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the crdinates and gives back the
            information at this specific position.
        """
        res = self._interface.get(request)
        # in the case of ab initio/DB add the masses
        if 'mass' in request and self.mode == 'ab initio':
            res['mass'] = self.get_masses()
        if 'atoms' in request and self.mode == 'ab initio':
            res['atoms'] = self.refgeo['atoms']
        return res

    def _get_refgeo(self, filename):
        natoms, atoms, crds = read_geom(refgeo_path)
        return natoms, {'atoms': atoms, 'crd': np.array(crds)}

    def get_masses(self):
        masses = []
        for i in range(self.natoms):
            # masses are given in an array of shape (natoms, 3) like
            # the coordinates so that they can be easily used in the
            # surface hopping algorithm
            masses += [[atomic_masses[atomname_to_id[self.refgeo['atoms'][i]]],
                        atomic_masses[atomname_to_id[self.refgeo['atoms'][i]]],
                        atomic_masses[atomname_to_id[self.refgeo['atoms'][i]]]]]
        return np.array(masses)
