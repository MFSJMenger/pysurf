import copy

import numpy as np

from . import engine

@engine.register_action
def read_npfile(filename:"file") -> "array2D":
    return np.loadtxt(filename)

@engine.register_action
def array2D_extract_columns(data: "array2D", columns: "list") -> "array2D":
    newshape = list(data.shape)
    newshape[1] = len(columns)
    datanew = np.empty(newshape)
    for idx, c in enumerate(columns):
        datanew[:,idx] = copy.deepcopy(data[:,c])
    return datanew

@engine.register_action
def array2D_extract_column(data: "array2D", column: "int") -> "array1D":
    datanew = copy.deepcopy(data[:, column])
    return datanew

@engine.register_action
def array_extract_column(data: "array", column: "int") -> "array1D":
    datanew = copy.deepcopy(data[:, column])
    return datanew

@engine.register_action
def array2D_add_column(data: "array2D", add: "array1D", append: "bool"=True) -> "array2D":
    append = False
    newshape = list(data.shape)
    newshape[1] += 1
    datanew = np.empty(newshape)
    if append is True:
        datanew[:, :-1] = data
        datanew[:, -1] = add
    else:
        datanew[:, 0] = add
        datanew[:, 1:] = data
    return datanew

@engine.register_action
def array2D_sort(data: "array2D") -> "array2D":
    data = data[np.argsort(data[:, 0])]
    return data

@engine.register_action
def array2D_add_ref(data: "array2D", columns: "list", ref: "float"):
    if isinstance(columns, list):
        for c in columns:
            data[:,c] = data[:,c] + ref
    if isinstance(columns, int):
        data[:, columns] = data[:, columns] + ref

@engine.register_action
def array2D_scale(data: "array2D", columns: "list", scale: "float"):
    for c in columns:
        data[:,c] = data[:,c] * scale
