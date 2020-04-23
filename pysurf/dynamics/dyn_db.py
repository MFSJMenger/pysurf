import numpy as np

from pysurf.database import Database
from pysurf.database.dbtools import DBVariable

class DynDB(Database):

    @staticmethod
    def _settings(natoms, nstates, model=False, nmodes=0):

        variables = {}
        variables['curr_state'] = DBVariable(int, ('frame', 'one'))
        variables['energy'] = DBVariable(np.double, ('frame', 'nstates'))
        variables['ekin'] = DBVariable(np.double, ('frame', 'one'))
        variables['epot'] = DBVariable(np.double, ('frame', 'one'))
        variables['etot'] = DBVariable(np.double, ('frame', 'one'))
        #check whether model or abinit calculations
        if model == False:
            variables['crd'] = DBVariable(np.double, ('frame', 'natoms', 'three'))
            variables['gradient'] = DBVariable(np.double, ('frame', 'nstates', 'natoms', 'three'))
            variables['mass'] = DBVariable(np.double, ('natoms'))
            variables['veloc'] = DBVariable(np.double, ('frame','natoms', 'three'))
            # more things can be stored if necessary
            dct = {'dimensions':{'frame': 'unlimited',
                                 'natoms': natoms,
                                 'nstates': nstates,
                                 'three': 3,
                                 'one': 1},
                   'variables': variables
                   }

        else:
            variables['crd'] = DBVariable(np.double, ('frame', 'nmodes'))
            variables['veloc'] = DBVariable(np.double, ('frame', 'nmodes'))
            variables['mass'] = DBVariable(np.double, ('nmodes',))
            variables['gradient'] = DBVariable(np.double, ('frame', 'nstates', 'nmodes'))
            # more things can be stored if necessary
            dct = {'dimensions':{'frame': 'unlimited',
                                 'nmodes': nmodes,
                                 'nstates': nstates,
                                 'three': 3,
                                 'one': 1},
                   'variables': variables
        }

        return dct    
    
    @classmethod
    def append(db, data, veloc, curr_state, ekin, epot, etot):
        db.append('crd', data['crd'])
        db.append('gradient', data['gradient'])
        db.append('energy', data['energy'])
        db.append('veloc', veloc)
        db.append('curr_state', curr_state)
        db.append('ekin', ekin)
        db.append('epot', epot)
        db.append('etot', etot)
        db.increase

        with open('prop.out','a') as outfile:
            outfile.write('{0} {1} {2}\n'.format(ekin, epot, etot))
            
    @classmethod
    def from_dynamics(cls, dbfile):
        return cls.load_db(dbfile)

    @classmethod
    def create_db(cls, dbfile, natoms, nstates, model, nmodes):
        return cls(dbfile, cls._settings(natoms, nstates, model, nmodes))


    def add_mass(self, mass):
        self.set('mass', mass)

    @property
    def len(self):
        return len(self['crd'])
