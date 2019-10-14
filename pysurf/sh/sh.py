import numpy as np
from ..random.shrandom import RandomNumberGeneratorNP


random = RandomNumberGeneratorNP()

def vv_xstep(crd, v, a, dt):
    crd += v*dt + 0.5*a*dt*dt
    return crd


def vv_vstep(v, a_old, a, dt):
    v += 0.5 * (a_old + a) * dt
    return v


def get_acceleration(g, m):
    return -g/m


def abs_gt(a, b):
    if  isinstance(a, float) or isinstance(a, int):
        return abs(a) > abs(b)
    return any(a[i] > b[i] for i in range(len(a)))


def compute_landau_zener_probability(dt, d_curr,
                                     d_prev_step, d_two_step_prev):
    """compute landau zener sh probability between two states"""
    # compute the second derivative of the energy difference by time
    finite_difference_grad = ((d_curr + d_two_step_prev - 2 * d_prev_step)
                              / dt**2)
    # compute the hopping probability
    return np.exp((-np.pi/2.0)  # pi/2
                  * np.sqrt(d_prev_step**3
                            / finite_difference_grad))


def compute_diff(e, i, j):
    """return e[i] - e[j]"""
    return e[i] - e[j]


def landau_zener_hop(prop):
    """Decide if Landom Zener Hop appears"""
    if (prop > random.get()):
        return True
    return False


class LandauZener(object):

    @staticmethod
    def landau_zener_select(iactive, nstates, dt,
                            e_curr, e_prev_step, e_two_step_prev):
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
            d_two_step_prev = compute_diff(e_two_step_prev, iactive, istate)
            d_prev_step = compute_diff(e_prev_step, iactive, istate)
            d_curr = compute_diff(e_curr, iactive, istate)
            # compute probability
            if (abs_gt(d_curr, d_prev_step) and
                    abs_gt(d_two_step_prev, d_prev_step)):
                curr_hop_prob = compute_landau_zener_probability(
                                    dt, d_curr, d_prev_step, d_two_step_prev)
                if curr_hop_prob > prop:
                    prop = curr_hop_prob
                    iselected = istate
                    print(f"{iselected} = {istate}")
        #
        if iselected is None:
            return None
        if landau_zener_hop(prop):
            return iselected
        return None

    @classmethod
    def landau_zener(cls, iactive, nstates, dt, e_curr,
                     e_prev_step, e_two_step_prev):

        return cls.landau_zener_select(iactive, nstates, dt, e_curr,
                                       e_prev_step, e_two_step_prev)


def landau_zener_surfacehopping(iactive, nsteps, random_seed, inp):

    random = RandomNumberGenerator(random_seed)

    for istep in range(nsteps):
        # 1) write step info
        # 2) call interface
        crd = vv_xstep(crd, v, a, dt)
        e, g = interface.get_qm(crd)
        # update acceleration
        a_old = a
        a = get_acceleration(g, m)
        # 
        v = vv_vstep(v, a_old, a, dt)
        # Landau Zener:
        iselected = LandauZener.landau_zener(iactive, nstates, dt, 
                e_curr, e_prev_step, e_two_step_prev)
        if iselected is not None:
            # change active state
            iactive = iselected
            # compute energy and gradient
            e, g = interface.get_qm(crd)
            # get acceleration
            a = get_acceleration(g, m)
