from pysurf.spp import SurfacePointProvider
from pysurf.database import PySurfDB

from . import engine

@engine.register_action
def crds_model(sampling:"spp", mode:"int", start:"float", stop:"float", npoints"int")
    spp.

