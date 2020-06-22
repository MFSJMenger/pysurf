from . import engine

engine.add_subtypes("anything", "request")
engine.add_subtypes("anything", "str")
engine.add_subtypes("anything", "crds")
engine.add_subtypes("anything", "crd")
engine.add_subtypes("anything", "array")

engine.add_subtypes("file", "str")

engine.add_subtypes("crd", "list")

engine.add_subtypes("array", "array2D")
engine.add_subtypes("array", "array1D")
engine.add_subtypes("array", "crds")
engine.add_subtypes("array", "crd")

engine.add_subtypes("list", "ilist")
