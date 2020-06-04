"""
PySurf Module:
    Validation and Training of Interpolators

Provide infrastructure for the training of interpolators
and test them against a validation set
"""
import numpy as np

from pysurf.database import PySurfDB
from pysurf.spp import SurfacePointProvider
from pysurf.logger import get_logger
from pysurf.colt import Colt

from scipy.optimize import minimize


class Validation(Colt):

    _questions = """
    [validate]
    db =
    properties = :: list
    optimize = False :: bool
    """

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_block("training", Training.questions)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        self.inter = Training.from_config(config['training'])
        if config['validate']['optimize'] is False:
            self.inter.validate(config['validate']['db'], config['validate']['properties'])
        else:
            self.inter.optimize(config['validate']['db'], config['validate']['properties'])


class Training(Colt):

    _questions = """
        spp = spp.inp :: existing_file
        savefile = :: str
    """

    @classmethod
    def from_config(cls, config):
        return cls(config['spp'], config['savefile'])

    def __init__(self, sppinp, savefile=None):
        #
        self.logger = get_logger('validate.log', 'validation', [])
        #
        config = self._get_spp_config(sppinp)
        #
        natoms, self.nstates, properties = self._get_db_info(config['use_db']['database'])
        self.spp = SurfacePointProvider(None, properties, self.nstates, natoms, None,
                                        logger=self.logger, config=config)
        #
        self.interpolator = self.spp.interpolator
        #
        self.savefile = savefile
        self.interpolator.train(savefile)

    def _get_spp_config(self, filename):
        questions = SurfacePointProvider.generate_questions(presets="""
                use_db=yes :: yes
                [use_db(yes)]
                fit_only = yes :: yes
                write_only = no :: no
                """)
        return questions.ask(config=filename, raise_read_error=False)

    def _get_db_info(self, database):
        db = PySurfDB.load_database(database, read_only=True)
        rep = db.dbrep
        natoms = rep.dimensions.get('natoms', None)
        if natoms is None:
            natoms = rep.dimensions['nmodes']
        nstates = rep.dimensions['nstates']
        return natoms, nstates, db.saved_properties

    def validate(self, filename, properties):
        db = PySurfDB.load_database(filename, read_only=True)
        self._compute(db, properties)

    def _compute(self, db, properties):
        norm = {prop: [] for prop in properties}
        ndata = len(db)

        for i, crd in enumerate(db['crd']):
            result = self.spp.request(crd, properties)
            #
            for prop in properties:
                norm[prop].append(np.copy(result[prop] - db[prop][i]))

        for name, value in norm.items():
            errors = self.compute_errors(name, value, ndata)
        return errors

    def compute_errors(self, name, prop, nele):
        prop = np.array(prop)
        #
        mse = np.sum(prop)/nele
        mae = np.sum(np.absolute(prop))/nele
        rmsd = np.sqrt(np.sum(prop**2)/nele)
        #
        maxval = np.amax(prop)
        minval = np.amin(prop)
        self.logger.info(f"{name}:\n mse = {mse}\n mae = {mae}\n"
                         f" rmsd = {rmsd}\n maxval = {maxval}\n minval={minval}\n")
        return {'mse': mse, 'mae': mae, 'rmsd': rmsd, 'max_error': maxval}

    def optimize(self, filename, properties):
        db = PySurfDB.load_database(filename, read_only=True)

        def _function(epsilon):
            self.interpolator.epsilon = epsilon
            self.interpolator.train()
            error = self._compute(db, properties)
            return error['rmsd']
        res = minimize(_function, self.interpolator.epsilon, options={
                       'xatol': 1e-8, 'disp': True})
        print(res)
        self.interpolator.epsilon = res.x[0]
        self.interpolator.train(self.savefile)


def eucl_norm(x, y):
    return np.linalg.norm(x - y)


if __name__ == '__main__':
    Validation.from_commandline()
