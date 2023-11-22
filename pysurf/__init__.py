# -*- coding: utf-8 -*-

"""Top-level package for pysurf."""

__author__ = """Maximilian F.S.J. Menger, Johannes Ehrmaier"""
__email__ = 'maximilian.menger@pci.uni-heidelberg.de'
__version__ = '0.1.0'

#
import os
#
from colt import PluginLoader
from .spp.spp import SurfacePointProvider, get_spp
from .spp import  AbinitioBase, Model, Interpolator
# load plugins
base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
core_plugins = os.path.join(base, "core_plugins")
user_plugins = os.path.join(base, "plugins")
# load core plugins
PluginLoader(core_plugins, ignorefile='plugins.ini')
# load user plugins
PluginLoader(user_plugins, ignorefile='plugins.ini')
