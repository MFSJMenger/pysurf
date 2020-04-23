import numpy as np

from pysurf.spp import SurfacePointProvider
#from pysurf.dynamics import *
from pysurf.database import Database
from pysurf.database.dbtools import DBVariable

from .base_propagator import PropagatorBase






class LandauZener(PropagatorBase):
    
    properties = ['energy', 'gradient'] 
    
    def run(self, nsteps, dt):
        e_curr = None
        e_prev_step = None
        e_two_step_prev = None
        
        if self.restart is False:
            iactive = int(self.init.state)
            print('Johannes condition:', self.init)
    
            # set random number
            self.init_random()
            # set starting crds
            crd = self.init.crd
            #
            v = self.init.veloc
            # start
            print('Johannes first crd: ', crd)
            data = self.get_data(crd, iactive)
            print('Johannes data:', data)
            nstates = len(data['energy'])
            print('Johannes: ', iactive)
            print('Johannes: ', data['gradient'][iactive], self.mass)
            a = get_acceleration(data['gradient'][iactive], self.mass)
            print('Johannes gradient: ', data['gradient'][iactive])
            print('Johannes all: ', crd, v, a, dt)
            #
            e_curr = data['energy']
            #
            for istep in range(2):
                # If not restart, first 2 steps are just to save energy!
                crd = vv_xstep(crd, v, a, dt)
                print('Johannes before get_data:', crd)
                data = self.get_data(crd, iactive) 
                print('Johannes after get_data:', crd)
                # update acceleration
                a_old = a
                a = get_acceleration(data['gradient'][iactive], self.mass)
                v = vv_vstep(v, a_old, a, dt)
                if e_prev_step is None:
                    e_prev_step = e_curr
                    e_curr = data['energy']
                    continue
                e_two_step_prev = e_prev_step
                e_prev_step = e_curr
                e_curr = data['energy']
                
            # check whether ab initio calculation
            eref = 0
    
        for istep in range(nsteps):
            # 1) write step info
            # 2) call interface
            crd = vv_xstep(crd, v, a, dt)
            data = self.get_data(crd, iactive)
            #
            # update acceleration
            a_old = a
            a = get_acceleration(data['gradient'][iactive], self.mass)
            # 
            v = vv_vstep(v, a_old, a, dt)
            # Landau Zener:
            iselected = LandauZener.landau_zener(iactive, nstates, dt, 
                    e_curr, e_prev_step, e_two_step_prev)
            if iselected is not None:
                # change active state
                dE = data['energy'][iselected] - data['energy'][iactive]
                ekin = calc_ekin(self.mass, v)
                if (dE > ekin):
                    print(f"Too few energy -> no hop")
                else:
                    iactive = iselected
                    print(f"new selected state: {iselected}")
                    # rescale velocity
                    v = rescale_velocity(ekin, dE, v)
                    # get acceleration
                    data = self.get_data(crd, iactive)
                    a = get_acceleration(data['gradient'][iactive], self.mass)
            ekin = calc_ekin(self.mass, v)
            epot = e_curr[iactive]
            etot = ekin + epot
    
            # If total energy changes too much stop the dynamics
            # This is likely to happen for interpolated data
            if istep == 0:
                eref = etot
            else:
                if eref < 0.5*etot or eref > 2.*etot:
                    break
            print(iactive, ekin, epot, etot)
            self.db.append(data, v, iactive, ekin, epot, etot)

    def get_data(self, crd, gradstate):
        res = self.spp.request(crd, self.properties, states=[gradstate])
        return res

def calc_ekin(masses, veloc):
    ekin = 0.0
    for i, mass in enumerate(masses.flatten()):
        ekin += 0.5*mass*veloc.flatten()[i]*veloc.flatten()[i]
    return ekin

def x_update(crd, v, a, dt):
    """Currently we only support velocity verlet"""
    crd += v*dt + 0.5*a*dt*dt
    return crd

def v_update(v, a_old, a, dt):
    v += 0.5 * (a_old + a) * dt
    return v

def get_acceleration(g, m):
    g = np.array(g)
    m = np.array(m)
    if g.shape == m.shape:
        return -g/m
    else:
        m_new = np.repeat(m, 3).reshape((len(m),3))
        return -g/m_new

def rescale_velocity(ekin, dE, v):
    factor = (1. - dE/ekin)
    print(f"factor = {factor}")
    v *= np.sqrt(factor)
    return v

def _rescale_velocity_along_v(ekin, dE, v):
    factor = (1. - dE/ekin)
    v *= np.sqrt(factor)
    return v

def vv_xstep(crd, v, a, dt):
    crd += v*dt + 0.5*a*dt*dt
    return crd


def vv_vstep(v, a_old, a, dt):
    v += 0.5 * (a_old + a) * dt
    return v

