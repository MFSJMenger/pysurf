import numpy as np
#
from pysurf.database.dbtools import DBVariable
from .base_sampling import InitialCondition
from .sampling import SamplingBase


class InitialConditionsFactory(SamplingBase):
    """Factory to store intial condition samplers"""
    # setup questions
    _questions = 'inherited'
    subquestions: 'inherited'
    # setup plugin
    _plugins_storage = '_methods'
    _is_plugin_specialisation = True
    _is_plugin_factory = True
    _register_plugin = False
    #
    condition = InitialCondition

    @classmethod
    def _generate_subquestions(cls, questions):
        """ This class will not be inherited """
        # Update _questions from Sampling by adding an additional question
        questions.add_questions_to_block("""
            # State on which trajectories start
            initial state = 0 :: int
            """)
        questions.generate_cases("method", {name: method.questions
                                 for name, method in cls._methods.items()})


class InitialConditionsBase(InitialConditionsFactory):
    """Base Class for Initial Conditions Sampler"""

    def get_condition(self, idx):
        if idx >= self.nconditions:
            return None
        crd = self._db.get('crd', idx)
        veloc = self._db.get('veloc', idx)
        #state = self._db.get('state', idx)
        return InitialCondition(crd, veloc, None)

    @property
    def _settings(self):
        settings = super()._settings

        # Add initial state
        settings['variables']['state'] = DBVariable(np.double, ('frame', 'one'))

        # Add velocities
        if self.model is False:
            settings['variables']['veloc'] = DBVariable(np.double, ('frame', 'natoms', 'three'))
        if self.model is True:
            settings['variables']['veloc'] = DBVariable(np.double, ('frame', 'nmodes'))

        return settings

    @staticmethod
    def _write_condition(db, cond):
        db.append('crd', cond.crd)
        db.append('veloc', cond.veloc)
        db.append('state', cond.state)
        db.increase
