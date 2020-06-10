import copy

import numpy as np

from . import engine

@engine.register_action(["filename"], "np.array")
def nparray_read(filename):
    return np.loadtxt(filename)

@engine.register_action(["np.array", "columns"], "np.array")
def nparray_extract_columns(data, columns):
    if isinstance(columns, list):
        newshape = list(data.shape)
        newshape[1] = len(columns)
        datanew = np.empty(newshape)
        for idx, c in enumerate(columns):
            datanew[:,idx] = copy.deepcopy(data[:,c])
    if isinstance(columns, int):
        newshape = list(data.shape)
        newshape[1] = 1
        datanew = np.empty(newshape)
        datanew[:, 0] = copy.deepcopy(data[:, columns])
    return datanew


@engine.register_action(["np.array", "columns", "float"])
def nparray_add_ref(data, columns, ref):
    if isinstance(columns, list):
        for c in columns:
            data[:,c] = data[:,c] + ref
    if isinstance(columns, int):
        data[:, columns] = data[:, columns] + ref

@engine.register_action(["np.array", "columns", "float"])
def nparray_scale(data, columns, scale):
    if isinstance(columns, list):
        for c in columns:
            data[:,c] = data[:,c] * scale
    if isinstance(columns, int):
        data[:, columns] = data[:, columns] * scale
