import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#
from colt import Colt
from qctools.converter import energy_converter, length_converter, time_converter

from .plot import Plot

class Plot3D(Plot):
    """ Plotting class for pysurf results. It uses Matplotlib to plot and Colt to 
        handle the user input
    """
    _questions = 'inherited'
    extend_questions : 'inherited'

    @classmethod
    def _extend_questions(cls, questions):
        questions.add_questions_to_block("""
        z_units = au :: str :: [au, eV, fs, a.u., nm, ang, bohr, ps, ns]
        z_start =  :: float, optional
        z_end = :: float, optional
        z_label = :: str
        z_label_unit = True :: bool
    """)

    def __init__(self, config):
        super().__init__(config)
        if config['matplotlib style'] is None:
            #Hard coded matplotlib stylesheet
            style = os.path.dirname(os.path.abspath(__file__))
            style = os.path.join(style, 'pysurf3d.mplstyle')
        else:
            style = config['matplotlib style']
        try:
            plt.style.use(style)
        except:
            raise Exception(f"Matplotlib style {style} could not be set")



    def surface_plot(self, data, x_units_in, y_units_in, z_units_in, ax=None, line_props={}, save_plot=True, show_plot=True):
        """ 
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
        
        if z_units_in is not None:
            try:
                zconverter = globals()[z_units_in[0]+'_converter'].get_converter(z_units_in[1], config['z_units'])
            except:
                    raise InputException("Converter or units not known")
        else:
            zconverter = lambda x: x

        #Convert data
        convdatax = np.empty_like(data[0])
        convdatay = np.empty_like(data[1])
        convdataz = np.empty_like(data[2])
        convdatax = xconverter(data[0])
        convdatay = yconverter(data[1])
        convdataz = zconverter(data[2])

        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)#fig.add_subplot(111, projection='3d')

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

        lbl = r"{}".format(config['z_label'])
        if config['z_label_unit']:
            lbl += r" [{}]".format(config['z_units'])
        ax.set_zlabel(lbl)
        
        scale = False
        lim = [config['z_start']]
        lim += [config['z_end']]
        if None in lim: scale = True
        ax.set_zlim(*lim, auto=scale)

        if config['z_units'] == 'a.u.':
           ax.set_zticklabels('' for z in ax.get_zticks())

        #plot data
        ax.plot_surface(convdatax, convdatay, convdataz, antialiased=True, shade=True, linewidth=0)

        #legend
        if config['legend'] is not None:
            ax.legend(*config['legend'])

        #reduce number of ticks
        plt.locator_params(nbins=4)
        
        ax.dist = 13
        #save plot
        if config['save_plot'] and save_plot:
            plt.savefig(config['save_plot']['plot_file'])
        
        #show plot
        if config['show_plot'] and show_plot:
            plt.show()

        return ax

