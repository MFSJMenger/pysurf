import itertools
import copy

import numpy as np
import matplotlib.pyplot as plt

from pysurf.utils import exists_and_isfile
from pysurf.analysis import Plot

from . import engine


@engine.register_action(["str"], "plot")
def setup_lineplot(plot_input):
    if exists_and_isfile(plot_input):
       plot_config = Plot.generate_input(plot_input, config=plot_input)
    else:
       plot_config = Plot.generate_input(plot_input)
    return [Plot(plot_config), None]

@engine.register_action(["plot","np.array", "style"], "plot")
def add_plot(plot, data, style):
    ax = plot[1]
    save_plot = False
    for i, style in zip(range(1, len(data[0])), itertools.cycle(style)):
        if i == (len(data[0])-1): save_plot = True
        ax = plot[0].line_plot(data[:,[0,i]], x_units_in=None, y_units_in=None, ax=ax, show_plot=False, save_plot=save_plot, line_props=style)
    plot[1] = ax
    return [plot[0], ax]

"""
Type style:
    style is a a list of dictionaries that can be given to the plot functions.
    The form of the dictionary is not further specified.
"""
@engine.register_action(["style", "style"], "style")
def combine_plotstyles(style1, style2):
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

@engine.register_action(output_typ="style")
def colors_standard():
    return [{'color': 'black'},
            {'color': [253/255, 86/255, 33/255]},
            {'color': [63/255, 81/255, 181/255]},
            {'color': [138/255,193/255,73/255]},
            {'color': [255/255, 233/255, 75/255]}, 
            {'color': [157/255, 157/255,157/255]}]

@engine.register_action(output_typ="style")
def linestyles_standard():
    return [{'linestyle': 'solid'},
            {'linestyle': 'dashed'},
            {'linestyle': 'dotted'}]

@engine.register_action(output_typ="style")
def linestyle_dotted():
    return [{'linestyle': 'dotted'}]

@engine.register_action(output_typ="style")
def linestyle_dashed():
    return [{'linestyle': 'dashed'}]

@engine.register_action(output_typ="style")
def linestyle_solid():
    return [{'linestyle': 'solid'}]

@engine.register_action(["plot"], None)
def show_plot(plot):
    plot[0] = plt.gcf()
    plt.show()

@engine.register_action(["plot", "file"], None)
def save_plot(plot, filename):
    plot[0].savefig(filename)
