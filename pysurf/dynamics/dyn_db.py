import numpy as np

from pysurf.database import PySurfDB

class DynDB(PySurfDB):
    variables = ['crd_equi', 'modes_equi', 'model', 'atomids', 'freqs_equi', 'masses', 'currstate', 'crd', 'veloc', 'energy', 'ekin', 'epot', 'etot']

    @classmethod
    def from_dynamics(cls, dbfile):
        info = info_database(dbfile)
        return cls.load_database(dbfile, info['variables'], info['dimensions'])

    @classmethod
    def create_db(cls, dbfile, sampling, nstates, props):
        variables = cls.variables
        for prop in props:
            if prop not in variables:
                variables += [prop]
        dims = sampling.info['dimensions']
        dims['nstates'] = nstates
        db = cls.generate_database(dbfile, variables, dims, model=sampling.model, sp=False)
        db.add_reference_entry(sampling.molecule, sampling.modes, sampling.model)
        return db

    def add_step(self, data, veloc, currstate, ekin, epot, etot):
        for entry in data:
            print('Johannes:', entry)
            print('Johannes:', data[entry])
            if entry == 'gradient':
                for state in entry:
                    self.append(entry, state, entry[state])
                continue
            self.append(entry, data[entry])
        self.append('currstate', currstate)
        self.append('ekin', ekin)
        self.append('epot', epot)
        self.append('etot', etot)
        self.increase


