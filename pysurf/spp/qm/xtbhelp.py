import numpy as np
#
from qctools import generate_filereader
from qctools import Event
from qctools import register_event_type


def end_of_loop():
    raise StopIteration


def read_xtb_gradient(iterator, NAtoms=0):
    
    # skip first NAtoms+1 lines
    for _ in range(NAtoms+2):
        next(iterator)
    # read next NAtoms lines and map it to a np.array
    out = np.array(list(list(map(float, next(iterator).replace('D','E').split())) 
                   for _ in range(NAtoms))).flatten(), 1
    print(out[0])
    return out
# add new event xtb_gradient
register_event_type('xtb_gradient', [{'NAtoms': int}, [], read_xtb_gradient])


def read_xtb_hessian(iterator, NAtoms=0):
    #skip first line
    next(iterator)
    # return np array of result
    return np.array(list(list(map(float, line.split())) for line in iterator)).flatten(), 1

# add new event xtb_hessian
register_event_type('xtb_hessian', [{'NAtoms': int}, [], read_xtb_hessian])

# define Events
Dipole = Event('Dipole', 
        'grep', {'keyword': 'molecular dipole:',
                 'ilen': 1,
                 'ishift': 3},
        func='split',
        func_kwargs={'idx': [1, 2, 3], 
                     'typ': [float, float, float]}
)
# 
Energy = Event('Energy', 
        'grep', {'keyword': 'TOTAL ENERGY',
                 'ilen': 1,
                 'ishift': 0},
        func='split',
        func_kwargs={'idx': 3, 'typ': float}
)
#
Gradient = Event('Gradient', 'xtb_gradient', {'NAtoms': 'NAtoms'})
#
Hessian = Event('Hessian', 'xtb_hessian', {'NAtoms': 'NAtoms'})

xtb_config = {
        'Energy': Energy,
        'Dipole': Dipole,
        'Gradient': Gradient,
        'Hessian': Hessian
}

XTBReader = generate_filereader('XTBReader', xtb_config)
