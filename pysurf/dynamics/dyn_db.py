import numpy as np

from pysurf.database import PySurfDB

class DynDB(PySurfDB):
    variables_molecule = ['crd_equi', 'modes_equi', 'model', 'atomids', 'freqs_equi', 'masses', 'currstate', 'crd', 'veloc', 'energy', 'ekin', 'epot', 'etot', 'time']
    variables_model =  ['crd_equi', 'modes_equi', 'model', 'freqs_equi', 'masses', 'currstate', 'crd', 'veloc', 'energy', 'ekin', 'epot', 'etot', 'time']

    @classmethod
    def from_dynamics(cls, dbfile):
        info = cls.info_database(dbfile)
        if 'atomids' in info['variables']:
            model = False
        else:
            model = True
        return cls.load_database(dbfile, info['variables'], info['dimensions'], model=model)

    @classmethod
    def create_db(cls, dbfile, sampling, nstates, props):
        if sampling.model:
            variables = cls.variables_model
        else:
            variables = cls.variables_molecule

        # Add additionally requested properties    
        for prop in props:
            if prop not in variables:
                variables += [prop]
        dims = sampling.info['dimensions']
        dims['nstates'] = nstates
        dims['nactive'] = 1
        db = cls.generate_database(dbfile, variables, dims, model=sampling.model, sp=False)
        db.add_reference_entry(sampling.system, sampling.modes, sampling.model)
        return db

    def add_step(self, time, data, veloc, currstate, ekin, epot, etot):
        self.append('time', time)
        for key, value in data.iter_data():
            if key == 'gradient':
                self.append(key, value[currstate])
                continue
            self.append(key, value)
        self.append('crd', data.crd)
        self.append('currstate', currstate)
        self.append('ekin', ekin)
        self.append('epot', epot)
        self.append('etot', etot)
        self.append('veloc', veloc)
        self.increase


