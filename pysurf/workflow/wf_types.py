from . import engine

engine.add_subtypes("file", "str")

engine.add_subtypes("crd", "list")

engine.add_subtypes("array", "array2D")
engine.add_subtypes("array", "array1D")
engine.add_subtypes("array", "crds")
engine.add_subtypes("array", "crd")
