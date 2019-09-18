import numpy as np


def vv_xstep(crd, v, a, dt):
    crd += v*dt + 0.5*a*dt*dt
    return crd


def vv_vstep(v, a_old, a, dt):
    v += 0.5 * (a_old + a) * dt
    return v


def surface_hopping_step(crd, v, a):
    pass

