import numpy as np

from pysurf.analysis import Plot
from pysurf.dynamics import DynDB
from pysurf.utils import exists_and_isfile
#
from colt import Colt

class PESTraj(Colt):
    _questions = """
    prop_db = prop.db :: existing_file
    
    #reference energy in atomic units
    reference_energy = :: float
    """

    _plot_input = 'plot_pes_traj.inp'

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        db = DynDB.load_database(config['prop_db'], read_only=True)
        nstates = db.nstates
        npoints = len(db)
        data = np.empty((npoints, nstates+1), dtype=float)
        data[:, 0] = np.array(db['time'])[:, 0]
        data[:, 1:] = np.array(db['energy'])-config['reference_energy']
        currstate = np.array(db['currstate'])[:, 0]
        plot_config = {}
        plot_config['x_label'] = 'time'
        plot_config['y_label'] = 'energy'
        plot_config['y_units'] = 'eV'
        plot_config['x_units'] = 'fs'
        plot_config['y_label_unit'] = True
        plot_config['x_label_unit'] = True
        plot_config = {'': plot_config}
#        plot_config['save_plot'] = {'plot_file': 'pes_traj.png'}

        presets = """
        x_label = time 
        y_label = energy
        y_units = eV
        x_units = fs
        y_label_unit = True
        x_label_unit = True
        """

        if exists_and_isfile(self._plot_input):
            plot_config = Plot.generate_input(self._plot_input, config=self._plot_input, presets=presets)
        else:
            plot_config = Plot.generate_input(self._plot_input, config=plot_config, presets=presets)

        myplot = Plot(plot_config)


        if nstates == 1:
            save_plot = True
            show_plot = True
        else:
            save_plot = False
            show_plot = False
        myax = myplot.line_plot(data[:,[0,1]], x_units_in=['time', 'au'], y_units_in=['energy', 'au'], ax=None, save_plot=save_plot, show_plot=show_plot)
        for state in range(1, nstates):
            myax = myplot.line_plot(data[:,[0, state+1]],  x_units_in=['time', 'au'], y_units_in=['energy', 'au'], ax=myax, show_plot=False, save_plot=False)
        
        curr_plot = np.empty((npoints, 2), dtype=float)
        curr_plot[:, 0] = data[:, 0]
        for idx, state in enumerate(currstate):
            curr_plot[idx, 1] = data[idx, int(state + 1)]
        myax = myplot.line_plot(curr_plot,  x_units_in=['time', 'au'], y_units_in=['energy', 'au'], ax=myax, save_plot=True, show_plot=True, line_props={'marker': 'o', 'color': 'red'})

if __name__ == "__main__":
    PESTraj.from_commandline()
