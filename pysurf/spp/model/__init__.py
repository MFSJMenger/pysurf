from ..methodbase import ModelBase
from .harmonic_oscillator import HarmonicOscillator
from .pyrazine_schneider import PyrazineSchneider
from .pyrazine_sala import PyrazineSala

plugins = [HarmonicOscillator, PyrazineSchneider, PyrazineSala]
