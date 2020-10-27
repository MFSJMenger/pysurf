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
from colt import Colt

from scipy.optimize import minimize


class Validation(Colt):

    _questions = """
    db =
    properties = :: list
    save_pes = __NONE__ :: str
    save_graddiff = __NONE__ :: str
    optimize = False :: bool
    """

    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_block("training", Training.questions)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        self.inter = Training.from_config(config['training'])
        #
        if config['optimize'] is False:
            self.inter.validate(config['db'], config['properties'])
        else:
            self.inter.optimize(config['db'], config['properties'])
        #
        if config['save_pes'] != '__NONE__':
            self.inter.save_pes(config['save_pes'], config['db'])

        if config['save_graddiff'] != '__NONE__':
            self.inter.save_graddiff(config['save_graddiff'], config['db'])

class Training(Colt):

    _questions = """
        spp = spp.inp :: existing_file
    """

    @classmethod
    def from_config(cls, config):
        return cls(config['spp'])

    def __init__(self, sppinp):
        #
        self.logger = get_logger('validate.log', 'validation', [])
        #
        config = self._get_spp_config(sppinp)
        #
        natoms, self.nstates, properties = self._get_db_info(config['use_db']['database'])
        atomids = [1 for _ in range(natoms)]
        self.spp = SurfacePointProvider.from_config(config, properties, self.nstates, natoms, 
                                                    atomids=atomids, logger=self.logger)
        #
        self.interpolator = self.spp.interpolator
        #
        self.weightsfile = self.spp.interpolator.weightsfile
        self.interpolator.train(self.weightsfile)

    def _get_spp_config(self, filename):
        questions = SurfacePointProvider.generate_questions(presets="""
                use_db=yes :: yes
                [use_db(yes)]
                write_only = no :: no
                [use_db(yes)::write_only(no)]
                fit_only = yes :: yes
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

    def save_graddiff(self, filename, database):
        db = PySurfDB.load_database(database, read_only=True)
        results, _ = self._compute(db, ['gradient'])

        def str_join(values):
            return ' '.join(str(val) for val in values)

        with open(filename, 'w') as f:
            graddiff = [np.sqrt(np.mean((fitted-exact)**2)) for (fitted, exact) in results['gradient']]
            f.write("\n".join(f"{i} {diff}" 
                              for i, diff in enumerate(graddiff)))

    def save_pes(self, filename, database):
        db = PySurfDB.load_database(database, read_only=True)
        results, _ = self._compute(db, ['energy'])

        def str_join(values):
            return ' '.join(str(val) for val in values)

        with open(filename, 'w') as f:
            f.write("\n".join(f"{i} {str_join(fitted)} {str_join(exact)}" 
                              for i, (fitted, exact) in enumerate(results['energy'])))

    def _compute(self, db, properties):
        norm = {prop: [] for prop in properties}
        ndata = len(db)

        for i, crd in enumerate(db['crd']):
            result = self.spp.request(crd, properties)
            #
            for prop in properties:
                if prop != 'gradient':
                    norm[prop].append([np.copy(result[prop]), np.copy(db[prop][i])])
                else:
                    norm[prop].append([np.copy(result[prop].data), np.copy(db[prop][i])])

        for name, value in norm.items():
            errors = self.compute_errors(name, value, ndata)

        return norm, errors

    def compute_errors(self, name, prop, nele):
        prop = np.array([val[0] - val[1] for val in prop])

        #
        mse = np.mean(prop)
        mae = np.mean(np.absolute(prop))
        rmsd = np.sqrt(np.mean(prop**2))
        rmsd_state = np.sqrt(np.mean(prop**2, axis=0))
        #
        maxval = np.amax(prop)
        minval = np.amin(prop)
        self.logger.info(f"{name}:\n mse = {mse}\n mae = {mae}\n"
                         f" rmsd = {rmsd}\n rmsd_state = {rmsd_state}\n maxval = {maxval}\n minval={minval}\n")
        return {'mse': mse, 'mae': mae, 'rmsd': rmsd, 'rmsd_state': rmsd_state, 'max_error': maxval}

    def optimize(self, filename, properties):
        db = PySurfDB.load_database(filename, read_only=True)

        def _function(epsilon):
            print('opt cycle', epsilon)
            self.interpolator.epsilon = epsilon[0]
            self.interpolator.train()
            _, error = self._compute(db, properties)
            print(error)
            return error['rmsd']

        res = minimize(_function, self.interpolator.epsilon, method='nelder-mead', tol=1e-4, options={
            'maxiter': 25, 'disp': True, 'xatol': 0.0001})
        print(res)
        self.interpolator.epsilon = res.x[0]
        self.interpolator.train(self.weightsfile)


def eucl_norm(x, y):
    return np.linalg.norm(x - y)


if __name__ == '__main__':
    Validation.from_commandline()
