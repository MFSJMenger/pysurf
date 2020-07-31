import itertools
import copy

import numpy as np
import matplotlib.pyplot as plt

from pysurf.utils import exists_and_isfile
from pysurf.utils.colors import COLORS
from pysurf.analysis import Plot

from . import engine


@engine.register_action
def setup_lineplot(plot_input: "file") -> "plot":
    if exists_and_isfile(plot_input):
       plot_config = Plot.generate_input(plot_input, config=plot_input)
    else:
       plot_config = Plot.generate_input(plot_input)
    return [Plot(plot_config), None]

@engine.register_action
def mesh_plot(plot_input: "file", data: "meshplotdata") -> "plot":
    if exists_and_isfile(plot_input):
       plot_config = Plot.generate_input(plot_input, config=plot_input)
    else:
       plot_config = Plot.generate_input(plot_input)

    plot = Plot(plot_config)
    ax = plot.mesh_plot(data, x_units_in=None, y_units_in=None, ax=None, show_plot=False, save_plot=True)
    return [plot, ax]

@engine.register_action
def add_plot(plot: "plot", data: "array2D", style: "style"=None):
    ax = plot[1]
    save_plot = False
    if style is not None:
        for i, style in zip(range(1, len(data[0])), itertools.cycle(style)):
            if i == (len(data[0])-1): save_plot = True
            ax = plot[0].line_plot(data[:,[0,i]], x_units_in=None, y_units_in=None, ax=ax, show_plot=False, save_plot=save_plot, line_props=style)
    else:
        for i in range(1, len(data[0])):
            if i == (len(data[0])-1): save_plot = True
            ax = plot[0].line_plot(data[:,[0,i]], x_units_in=None, y_units_in=None, ax=ax, show_plot=False, save_plot=save_plot)
    plot[1] = ax

"""
Type style:
    style is a a list of dictionaries that can be given to the plot functions.
    The form of the dictionary is not further specified.
"""
@engine.register_action
def combine_plotstyles(style1: "style", style2: "style") -> "style":
    """ This function combines two styles
    """
    nentries = max(len(style1), len(style2))
    style1 = itertools.cycle(style1)
    style2 = itertools.cycle(style2)
    style = []
    #combine the two dicts
    for i in range(nentries):
        dict1 = copy.deepcopy(next(style1))
        dict1.update(next(style2))
        style += [dict1]
    return style

@engine.register_action
def standard_colors() -> "style":
    return [{'color': COLORS['black']},
            {'color': COLORS['lightblue']},
            {'color': COLORS['orange']},
            {'color': COLORS['darkgrey']},
            {'color': COLORS['green']}, 
            {'color': COLORS['gold']}]

@engine.register_action
def linestyles_standard() -> "style":
    return [{'linestyle': 'solid'},
            {'linestyle': 'dashed'},
            {'linestyle': 'dotted'}]

@engine.register_action
def linestyle_dotted() -> "style":
    return [{'linestyle': 'dotted'}]

@engine.register_action
def linestyle_dashed() -> "style":
    return [{'linestyle': 'dashed'}]

@engine.register_action
def linestyle_solid() -> "style":
    return [{'linestyle': 'solid'}]

@engine.register_action
def show_plot(plot: "plot"):
    plot[0] = plt.gcf()
    plt.show()

@engine.register_action
def save_plot(plot: "plot", filename: "file"):
    plot[0].savefig(filename)
