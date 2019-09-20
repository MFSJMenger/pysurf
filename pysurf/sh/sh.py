import numpy as np


def vv_xstep(crd, v, a, dt):
    crd += v*dt + 0.5*a*dt*dt
    return crd


def vv_vstep(v, a_old, a, dt):
    v += 0.5 * (a_old + a) * dt
    return v


def abs_gt(a, b):
    return abs(a) > abs(b)

def compute_landau_zener_probability(dt, d_curr, 
        d_prev_step, d_two_step_prev):
    finite_difference_grad = (d_curr + d_two_step_prev - 2 * d_prev_step) / dt**2        
class LandauZener(object):

    @staticmethod
    def compute_diff(e, i, j):
        """return e[i] - e[j]"""
        return e[i] - e[j]

    @staticmethod
    def landau_zener_select(iactive, nstates, dt, 
            e_curr, e_prev_step, e_two_step_prev):
        """"takes in the (adiabatic?) energies at
            the current, the previous, and the two steps previous time step

            select which state is the most likely to hop to
        """

        for istate in range(nstates):
            # 
            if istate == iactive:
                continue
            # select state
            d_two_step_prev = compute_diff(e_two_step_prev, iactive, istate)
            d_prev_step = compute_diff(e_prev_step, iactive, istate)
            d_curr = compute_diff(e_curr, iactive, istate)
            if (abs_gt(d_curr, d_prev_step) and abs_gt(d_two_step_prev, d_prev_step)):
                chop = compute_landau_zener_probability(dt, d_curr, 
                        d_prev_step, d_two_step_prev)





