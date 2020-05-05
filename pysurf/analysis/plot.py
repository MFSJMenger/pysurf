import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from pysurf.colt import Colt
from pysurf.qctools.converter import energy_converter, length_converter, time_converter

class Plot(Colt):
    """ Plotting class for pysurf results. It uses Matplotlib to plot and Colt to 
        handle the user input
    """

    _questions = """
        x_units = au :: str :: [au, eV, fs, a.u., nm, ang, bohr, ps, ns]
        y_units = au :: str :: [au, eV, fs, a.u., nm, ang, bohr, ps, ns]
        x_start =  :: float, optional
        x_end = :: float, optional
        y_start =  :: float, optional
        y_end = :: float, optional
        y_label = :: str
        y_label_unit = True :: bool
        x_label = :: str
        x_label_unit = True :: bool
        legend = :: python(dict), optional
        title =  :: str, optional
        save_plot = yes :: str :: [yes, no]
        show_plot = True :: bool
        rcparams = :: python(dict), optional
        matplotlib style = :: str, optional
    """

    _save_plot = {
            'yes': "plot_file = plot.png :: file",
            'no': ""
            }

    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases('save_plot', {name: mode for name, mode in cls._save_plot.items()})

    @classmethod
    def from_inputfile(cls, inputfile):
        config = self.generate_input(inputfile, config=inputfile)
        return cls(config)

    def __init__(self, config):
        """
        Initialising the plotting class

        Args:
            config:
                Colt information according to the class questions
        """
        self.config = config
        if config['matplotlib style'] is None:
            #Hard coded matplotlib stylesheet
            style = os.path.dirname(os.path.abspath(__file__))
            style = os.path.join(style, 'pysurf.mplstyle')
        else:
            style = config['matplotlib style']
        try:
            plt.style.use(style)
        except:
            raise Exception(f"Matplotlib style {style} could not be set")

        if config['rcparams'] is not None:
            for key, value in config['rcparams'].items():
                mpl.rcParams[key] = value


    def line_plot(self, data, x_units_in, y_units_in, ax=None, line_props={}, save_plot=True, show_plot=True):
        """ 
        function that plots the data in a line plot

        Args:
            data, numpy array 
                size (:, 2), containing the x and y data for the plot.

            x_units_in, None or tuple:
                If it is None, x data will not be converted
                If it is a tuple of length 2, it has to contain the type of the x-data, i.e. length, energy or time
                and the units of the input data. The data will be transformed accordingly from this unit to the
                desired unit given in the Colt questins, i.e. self.config['x_units']

            y_units_in, None or tuple:
                If it is None, ydata will not be converted
                If it is a tuple of length 2, it has to contain the type of the y-data, i.e. length, energy or time
                and the units of the input data. The data will be transformed accordingly from this unit to the
                desired unit given in the Colt questions, i.e. self.config['y_units']

            ax, optional, pyplot axes object:
                If ax is set, the Line Plot will be added to the chart, otherwise a new figure is created for the plot.

            line_props, optional, dict:
                dictionary containing the properties like linestyl, color, width, ... for the line plot. It will be used as is.
                

        Output:
            ax, pyplot axes object:
                ax is the updated axes object with the new plot, or the newly created object if ax was not set.

        """


        config = self.config
         
        #Get converter for x-data
        if x_units_in is not None:
            try:
                xconverter = globals()[x_units_in[0]+'_converter'].get_converter(x_units_in[1], config['x_units'])
            except:
                raise InputException("Converter or units not known")
        else:
            xconverter = lambda x: x
        #Get converter for y-data
        if y_units_in is not None:
            try:
                yconverter = globals()[y_units_in[0]+'_converter'].get_converter(y_units_in[1], config['y_units'])
            except:
                    raise InputException("Converter or units not known")
        else:
            yconverter = lambda x: x
        
        #Convert data
        convdata = np.empty_like(data)
        convdata[:, 0] = xconverter(data[:, 0])
        convdata[:, 1] = yconverter(data[:, 1])
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

#        if config['title'] != None: ax.set_title(config['title'])


        #set axis and axis label
        lbl = r"{}".format(config['x_label'])
        if config['x_label_unit']:
            lbl += r" [{}]".format(config['x_units'])
        ax.set_xlabel(lbl)
 
        scale = False
        lim = [config['x_start']]
        lim += [config['x_end']]
        if None in lim: scale = True
        ax.set_xlim(*lim, auto=scale)

        if config['x_units'] == 'a.u.':
            ax.set_xticklabels('' for x in ax.get_xticks())
        
        lbl = r"{}".format(config['y_label'])
        if config['y_label_unit']:
            lbl += r" [{}]".format(config['y_units'])
        ax.set_ylabel(lbl)
        
        scale = False
        lim = [config['y_start']]
        lim += [config['y_end']]
        if None in lim: scale = True
        ax.set_ylim(*lim, auto=scale)

        if config['y_units'] == 'a.u.':
            ax.set_yticklabels('' for y in ax.get_yticks())

        #plot data
        ax.plot(convdata[:, 0], convdata[:, 1], **line_props)

        #legend
        if config['legend'] is not None:
            ax.legend(**config['legend'])

        #save plot
        if config['save_plot'] and save_plot:
            plt.savefig(config['save_plot']['plot_file'])

        #show plot
        if config['show_plot'] and show_plot:
            plt.show()

        return ax

