import numpy as np

from pysurf.spp.spp import SurfacePointProvider
from pysurf.sh.sh import *
from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable

#spp = SurfacePointProvider('./test.inp')
#out = spp.get({'coord': np.zeros(3, dtype=np.double), 'energy': None, 'gradient': None})

class Save(object):

    def __init__(self, filename, data, header=None):

        crds = data['coord']
            

        variables = {}
        variables['curr_state'] = DBVariable(int, ('frame', 'one'))
        variables['energy'] = DBVariable(np.double, ('frame', 'nstates'))
        variables['ekin'] = DBVariable(np.double, ('frame', 'one'))
        variables['epot'] = DBVariable(np.double, ('frame', 'one'))
        variables['etot'] = DBVariable(np.double, ('frame', 'one'))
        #check whether model or abinit calculations
        if len(crds.shape) == 2:

            variables['coord'] = DBVariable(np.double, ('frame', 'natoms', 'three'))
            variables['gradient'] = DBVariable(np.double, ('frame', 'nstates', 'natoms', 'three'))
            variables['mass'] = DBVariable(np.double, ('natoms','three'))
            variables['veloc'] = DBVariable(np.double, ('frame','natoms', 'three'))
            # more things can be stored if necessary
            dct = {'dimensions':{'frame': 'unlimited',
                                 'natoms': crds.shape[0],
                                 'nstates': len(data['energy']),
                                 'three': 3,
                                 'one': 1},
                   'variables': variables
                   }

        else:
            variables['coord'] = DBVariable(np.double, ('frame', 'nmodes'))
            variables['veloc'] = DBVariable(np.double, ('frame', 'nmodes'))
            variables['mass'] = DBVariable(np.double, ('nmodes',))
            variables['gradient'] = DBVariable(np.double, ('frame', 'nstates', 'nmodes'))
            # more things can be stored if necessary
            dct = {'dimensions':{'frame': 'unlimited',
                                 'nmodes': len(crds),
                                 'nstates': len(data['energy']),
                                 'three': 3,
                                 'one': 1},
                   'variables': variables
        }
        
        self.db = Database(filename, dct)

    def append(self, data, veloc, curr_state, ekin, epot, etot):
        self.db.append('coord', data['coord'])
        self.db.append('gradient', data['gradient'])
        self.db.append('energy', data['energy'])
        self.db.append('veloc', veloc)
        self.db.append('curr_state', curr_state)
        self.db.append('ekin', ekin)
        self.db.append('epot', epot)
        self.db.append('etot', etot)
        self.db.increase

        with open('prop.out','a') as outfile:
            outfile.write('{0} {1} {2}\n'.format(ekin, epot, etot))
            

    def add_mass(self, mass):
        self.db.set('mass', mass)

#    def __del__(self):
#        self.f.close()


def rescale_velocity(ekin, dE, v):
    factor = (1. - dE/ekin)
    print(f"factor = {factor}")
    v *= np.sqrt(factor)
    return v

def landau_zener_surfacehopping(init_cond, iactive, nsteps, random_seed, inp, dt, model=False):

#    mass = np.array([1.0, 1.0, 1.0])

#    save_crd = Save('crd.txt', "M1 M2 M3")
#    save_veloc = Save('veloc.txt', "M1 M2 M3")
#    save_acc = Save('acc.txt', "M1 M2 M3")
#    save_energy = Save("energy.txt", "iactive EKin EPot ETot")
#    save_pes = Save("pes.txt", "states")
    e_curr = None
    e_prev_step = None
    e_two_step_prev = None

    # init SPP
    spp = SurfacePointProvider(inp)
    # set random number
    random = RandomNumberGeneratorNP(random_seed)
    # set starting coords
    crd = init_cond.crd
    #
    v = init_cond.veloc
    # start
    data = get_data(spp, crd)
    nstates = len(data['energy'])
    a = get_acceleration(data['gradient'][iactive], data['mass'])
    #
    e_curr = data['energy']
    for istep in range(2):
        """If not restart, first 2 steps are just to save energy!"""
        crd = vv_xstep(crd, v, a, dt)
        data = get_data(spp, crd) 
        # update acceleration
        a_old = a
        a = get_acceleration(data['gradient'][iactive], data['mass'])
        v = vv_vstep(v, a_old, a, dt)
        if e_prev_step is None:
            e_prev_step = e_curr
            e_curr = data['energy']
            continue
        e_two_step_prev = e_prev_step
        e_prev_step = e_curr
        e_curr = data['energy']
        
    save = Save('prop.db', data)
    #check whether ab initio calculation
    save.add_mass(data['mass'])

    for istep in range(nsteps):
        # 1) write step info
        # 2) call interface
        crd = vv_xstep(crd, v, a, dt)
        data = get_data(spp, crd)
        #
        e_two_step_prev = e_prev_step
        e_prev_step = e_curr
        e_curr = data['energy']
        # update acceleration
        a_old = a
        a = get_acceleration(data['gradient'][iactive], data['mass'])
        # 
        v = vv_vstep(v, a_old, a, dt)
        # Landau Zener:
        iselected = LandauZener.landau_zener(iactive, nstates, dt, 
                e_curr, e_prev_step, e_two_step_prev)
        if iselected is not None:
            # change active state
            dE = data['energy'][iselected] - data['energy'][iactive]
            ekin = calc_ekin(data['mass'], v)
            if (dE > ekin):
                print(f"Too few energy -> no hop")
            else:
                iactive = iselected
                print(f"new selected state: {iselected}")
                # rescale velocity
                v = rescale_velocity(ekin, dE, v)
                # get acceleration
                a = get_acceleration(data['gradient'][iactive], data['mass'])
#        save_crd.save(" ".join("%12.8f" % c for c in crd.flatten()))
#        save_veloc.save(" ".join("%12.8f" % _v for _v in v.flatten()))
#        save_acc.save(" ".join("%12.8f" % _a for _a in a.flatten()))
        ekin = calc_ekin(data['mass'], v)
        epot = e_curr[iactive]
        etot = ekin + epot
        print(iactive, ekin, epot, etot)
        save.append(data, v, iactive, ekin, epot, etot)
#        save_energy.save("%8d %12.8f   %12.8f    %12.8f" % (iactive, ekin, epot, etot))
#        save_pes.save("%12.8f %12.8f   %12.8f" % (data['energy'][0], data['energy'][1], data['energy'][2]))

def calc_ekin(masses, veloc):
    ekin = 0.0
    for i, mass in enumerate(masses.flatten()):
        ekin += 0.5*mass*veloc.flatten()[i]*veloc.flatten()[i]
    return ekin

def get_data(spp, crd):
    res = spp.get({'coord': crd, 'mass': None, 'gradient': None, 'energy': None})
    return res
