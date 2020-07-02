# -*- coding: utf-8 -*-

"""Top-level package for pysurf."""

__author__ = """Maximilian F.S.J. Menger, Johannes Ehrmaier"""
__email__ = 'menger.maximilian@gmail.com'
__version__ = '0.1.0'

#
import os
#
from .colt import PluginLoader
from .spp.spp import SurfacePointProvider
from .spp import  AbinitioBase, Model, Interpolator
# load plugins
plugins = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../plugins")
PluginLoader(plugins)
