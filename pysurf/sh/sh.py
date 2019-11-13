import numpy as np
from .shbase import SurfaceHoppingBase


def get_data(spp, crd):
    res = spp.get({'crd': crd, 'mass': None, 'gradient': None, 'energy': None})
    return res


def calc_ekin(masses, veloc):
    ekin = 0.0
    for i, mass in enumerate(masses.flatten()):
        ekin += 0.5*mass*veloc.flatten()[i]*veloc.flatten()[i]
    return ekin


class LandauZener(SurfaceHoppingBase):

    def __init__(self, configfile):
        """setup surface hopping base"""
        super().__init__(self, configfile)

    def run(self):
        """actual surface hopping run"""
        # these here should got from init
        e_curr = self.storage['energy'].current
        e_prev_step =  self.storage['energy'].prev
        e_two_step_prev =  self.storage['energy'].twoprev
        # set starting crd
        crd = self.storage['crd']
        #
        v = self.storage['veloc']
        # start
        data = get_data(spp, crd)
        nstates = len(data['energy'])
        a = self.get_acceleration(data['gradient'][iactive], data['mass'])
        #
        if self.irestart is False:
            # setup
            e_curr = data['energy']
            for istep in range(2):
                # If not restart, first 2 steps are just to save energy!
                crd, v, a, data = self._propagation_step(crd, v, a, dt)
                if e_prev_step is None:
                    e_prev_step = e_curr
                    e_curr = data['energy']
                    continue
                e_two_step_prev = e_prev_step
                e_prev_step = e_curr
                e_curr = data['energy']
        else:
            raise Exception("not implemented!")
        # TODO: that should be removed!
        save = Save('prop.db', data)
        # check whether ab initio calculation
        save.add_mass(data['mass'])
        #
        for istep in range(nsteps):
            # 1) write step info
            # 2) call interface
            crd, v, a, data = self._propagation_step(crd, v, a, dt)
            # update energy
            e_two_step_prev = e_prev_step
            e_prev_step = e_curr
            e_curr = data['energy']
            # Landau Zener:
            iselected = self.landau_zener(iactive, nstates, dt, 
                    e_curr, e_prev_step, e_two_step_prev)
            if iselected is not None:
                # change active state
                dE = data['energy'][iselected] - data['energy'][iactive]
                ekin = calc_ekin(data['mass'], v)
                # TODO: all logging should be done automatically
                if (dE > ekin):
                    print(f"Too few energy -> no hop")
                else:
                    iactive = iselected
                    print(f"new selected state: {iselected}")
                    # rescale velocity
                    v = self._rescale_velocity_along_v(ekin, dE, v)
                    # get new acceleration
                    a = self.get_acceleration(data['gradient'][iactive], data['mass'])
            ekin = calc_ekin(data['mass'], v)
            epot = e_curr[iactive]
            etot = ekin + epot
            print(iactive, ekin, epot, etot)
            save.append(data, v, iactive, ekin, epot, etot)

    def _propagation_step(self, crd, v, a, dt):
        crd = self.x_update(crd, v, a, dt)
        data = get_data(spp, crd) 
        # update acceleration
        a_old = a
        a = self.get_acceleration(data['gradient'][iactive], data['mass'])
        v = self.v_update(v, a_old, a, dt)
        return crd, v, a, data

    def _init_trajectory(self):
        """setup a new trajectory run, returns if restart or not!"""
        self.init_random(10)
        return False

    def _restart_trajectory(self):
        raise Exception("not implemented yet")

    @classmethod
    def landau_zener_select(cls, iactive, nstates, dt, e_curr, e_prev_step, e_two_step_prev):
        """"takes in the (adiabatic?) energies at
            the current, the previous, and the two steps previous time step

            select which state is the most likely to hop to
        """
        iselected = None
        prop = -1.0
        for istate in range(nstates):
            #
            if istate == iactive:
                continue
            # compute energy differences
            d_two_step_prev = cls.compute_diff(e_two_step_prev, iactive, istate)
            d_prev_step = cls.compute_diff(e_prev_step, iactive, istate)
            d_curr = cls.compute_diff(e_curr, iactive, istate)
            # compute probability
            if (cls.abs_gt(d_curr, d_prev_step) and
                    cls.abs_gt(d_two_step_prev, d_prev_step)):
                curr_hop_prob = cls.compute_landau_zener_probability(
                                    dt, d_curr, d_prev_step, d_two_step_prev)
                if curr_hop_prob > prop:
                    prop = curr_hop_prob
                    iselected = istate
                    print(f"{iselected} = {istate}")
        #
        if iselected is None:
            return None
        if cls.landau_zener_hop(prop):
            return iselected
        return None

    @staticmethod
    def abs_gt(a, b):
        if  isinstance(a, float) or isinstance(a, int):
            return abs(a) > abs(b)
        return any(a[i] > b[i] for i in range(len(a)))

    @staticmethod
    def compute_landau_zener_probability(dt, d_curr, d_prev_step, d_two_step_prev):
        """compute landau zener sh probability between two states"""
        # compute the second derivative of the energy difference by time
        finite_difference_grad = ((d_curr + d_two_step_prev - 2 * d_prev_step)
                                  / dt**2)
        # compute the hopping probability
        return np.exp((-np.pi/2.0)  # pi/2
                      * np.sqrt(d_prev_step**3
                                / finite_difference_grad))

    @staticmethod
    def compute_diff(e, i, j):
        """return e[i] - e[j]"""
        return e[i] - e[j]

    def landau_zener_hop(self, prop):
        """Decide if Landom Zener Hop appears"""
        if (prop > self.random()):
            return True
        return False
