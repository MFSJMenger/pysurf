import numpy as np
from pysurf.random.shrandom import RandomNumberGeneratorNP

from pysurf.spp import SurfacePointProvider
#from pysurf.dynamics import *

from .base_propagator import PropagatorBase






class LandauZener(PropagatorBase):
    
    properties = ['energy', 'gradient'] 
    
    def run(self, nsteps, dt, seed=16661):
        self.dt = dt
        self.nsteps = nsteps
        self.e_curr = None
        self.e_prev_step = None
        self.e_two_step_prev = None
        
        # set random number
        self.init_random(seed)

        if self.restart is False:
            self.setup_new()
        else:
            self.setup_from_db()

        if self.start > self.nsteps:
            self.logger.info("Dynamics already performed for more steps!")
            return 

        for istep in range(self.start, nsteps+1):
            infotext = None
            # 1) write step info
            # 2) call interface
            self.crd = vv_xstep(self.crd, self.v, self.a, self.dt)
            data = self.call_spp()
            #
            # update acceleration
            self.a_old = self.a
            self.a = get_acceleration(data['gradient'][self.iactive], self.masses)
            # 
            self.v = vv_vstep(self.v, self.a_old, self.a, self.dt)

            # Landau Zener:
            self.iselected = self.lz_select()
            if self.iselected is not None:
                # change active state
                self.dE = data['energy'][iselected] - data['energy'][self.iactive]
                self.ekin = calc_ekin(self.masses, self.v)
                if (self.dE > self.ekin):
                    infotext = (f"*   Frustrated hop: Too few energy -> no hop\n")
                else:
                    infotext = (f"*   LandauZener hop: {self.iactive} -> {iselected}\n")
                    self.iactive = self.iselected
                    # rescale velocity
                    self.v = rescale_velocity(self.ekin, self.dE, self.v)
                    # get acceleration
                    data = self.call_spp()
                    self.a = get_acceleration(data['gradient'][self.iactive], self.masses)
            self.ekin = calc_ekin(self.masses, self.v)
            self.epot = self.e_curr[self.iactive]
            self.etot = self.ekin + self.epot
            
            time = dt * istep
            if istep == 0:
                dE = 0.
            else:
                dE = self.etot - etot_old
            self.log_step(istep, time, self.iactive, dE, self.ekin, self.epot, self.etot, infotext) 
            self.db.add_step(data, self.v, self.iactive, self.ekin, self.epot, self.etot)
            etot_old = self.etot

    def call_spp(self, crd=None, gradstate=None):
        if crd is None:
            crd = self.crd
        if gradstate is None:
            gradstate = self.iactive
        res = self.spp.request(crd, self.properties, states=[gradstate])
        return res

    def setup_new(self):
            self.iactive = int(self.init.state)
            # set starting crds
            self.crd = self.init.crd
            # set starting veloc
            self.v = self.init.veloc
            # start
            data = self.call_spp()
            nstates = len(data['energy'])
            self.a = get_acceleration(data['gradient'][self.iactive], self.masses)
            #
            self.e_curr = data['energy']
            #
            for istep in range(2):
                # If not restart, first 2 steps are just to save energy!
                self.crd = vv_xstep(self.crd, self.v, self.a, self.dt)
                data = self.call_spp() 
                # update acceleration
                self.a_old = self.a
                self.a = get_acceleration(data['gradient'][self.iactive], self.masses)
                self.v = vv_vstep(self.v, self.a_old, self.a, self.dt)
                if self.e_prev_step is None:
                    self.e_prev_step = self.e_curr
                    self.e_curr = data['energy']
                    continue
                self.e_two_step_prev = self.e_prev_step
                self.e_prev_step = self.e_curr
                self.e_curr = data['energy']
            self.start = 0

    def setup_from_db(self):
        #two previous steps are needed
        if self.db.len < 2:
            return self.setup_new()
        else:
            self.crd = self.db.get('crd', -1)
            self.iactive = int(self.db.get('currstate', -1))
            self.v = self.db.get('veloc', -1)
            grad = self.db.get('gradient', -1)
            self.a = get_acceleration(grad, self.masses)
            self.e_curr = self.db.get('energy',-1)
            self.e_prev_step = self.db.get('energy', -2)
            self.e_two_step_prev = self.db.get('energy', -3)
            self.start = self.db.len            

    def lz_select(self):
        """"takes in the (adiabatic?) energies at
            the current, the previous, and the two steps previous time step
            select which state is the most likely to hop to
        """
        iselected = None
        prob = -1.0
        for istate in range(self.nstates):
            #
            if istate == self.iactive:
                continue
            # compute energy differences
            d_two_step_prev = compute_diff(self.e_two_step_prev, self.iactive, istate)
            d_prev_step = compute_diff(self.e_prev_step, self.iactive, istate)
            d_curr = compute_diff(self.e_curr, self.iactive, istate)
            # compute probability
            if (abs_gt(d_curr, d_prev_step) and abs_gt(d_two_step_prev, d_prev_step)):
                curr_hop_prob = _lz_probability(self, self.dt, d_curr, d_prev_step, d_two_step_prev)
                if curr_hop_prob > prob:
                    prob = curr_hop_prob
                    iselected = istate
        #
        if iselected is None:
            return None
        if prob > self.random():
            return iselected
        return None


    def _lz_probability(self, dt, d_curr, d_prev_step, d_two_step_prev):
        """compute landau zener sh probability between two states"""
        # compute the second derivative of the energy difference by time
        finite_difference_grad = ((d_curr + d_two_step_prev - 2 * d_prev_step) / dt**2)
        # compute the hopping probability
        return np.exp((-np.pi/2.0) * np.sqrt(d_prev_step**3 / finite_difference_grad))


    def init_random(self, seed=16661):
        self._random = RandomNumberGeneratorNP(seed)

    def random(self):
        return self._random.get()


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

def compute_diff(e, i, j):
    """return e[i] - e[j]"""
    return e[i] - e[j]

def abs_gt(a, b):
    return abs(a) > abs(b)



