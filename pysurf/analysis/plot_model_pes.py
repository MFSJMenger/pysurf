import numpy as np

from pysurf.analysis import Plot
from pysurf.dynamics import DynDB
from pysurf.utils import exists_and_isfile
from pysurf.spp import ModelBase, Request

from colt import Colt


class PESModel(Colt):
    _questions = """
    model = :: str
    mode = 0 :: int
    start_plot = -5 :: float
    end_plot = 5 :: float
    """

    _plot_input = "plot_pes_model.inp"
    
    @classmethod
    def _extend_questions(cls, questions):
        questions.generate_cases("model", {name: model.questions for name, model in ModelBase._models.items()})

    @classmethod
    def from_inputfile(cls, inputfile):
        config = cls.generate_input(inputfile, config=inputfile)
        return cls(config)

    @classmethod
    def from_config(cls, config):
        return cls(config)

    def __init__(self, config):
        model = ModelBase._models[config['model'].value].from_config(config['model'])
        nstates = model.nstates
        crd_equi = model.crd
        npoints = 100
        data = np.empty((npoints, nstates+1), dtype=float)
        data[:, 0] = np.linspace(config['start_plot'], config['end_plot'], npoints)
        for idx, point in enumerate(data[:, 0]):
            crd = np.zeros(model.crd.size)
            crd[config['mode']] = point
            data[idx, 1:] = model.get(Request(crd, properties=['energy'], states=[i for i in range(nstates)]))['energy']
        print(data)
        plot_config = {}
        plot_config['x_label'] = f"mode config['mode']"
        plot_config['x_start'] = config['start_plot']
        plot_config['x_end'] = config['end_plot']
        plot_config['y_label'] = 'energy'
        plot_config['y_units'] = 'eV'
        plot_config['x_units'] = 'au'
        plot_config['y_label_unit'] = True
        plot_config['x_label_unit'] = True
        plot_config = {'': plot_config}

        presets = """
        x_label = mode 
        y_label = energy
        y_units = eV
        x_units = au
        y_label_unit = True
        x_label_unit = True
        """

        if exists_and_isfile(self._plot_input):
            plot_config = Plot.generate_input(self._plot_input, config=self._plot_input, presets=presets)
        else:
            plot_config = Plot.generate_input(self._plot_input, config=None, presets=None)

        myplot = Plot(plot_config)


        if nstates == 1:
            save_plot = True
            show_plot = True
        else:
            save_plot = False
            show_plot = False
        myax = myplot.line_plot(data[:,[0,1]], x_units_in=['energy', 'a.u.'], y_units_in=['energy', 'au'], ax=None, save_plot=save_plot, show_plot=show_plot)
        for state in range(1, nstates-1):
            myax = myplot.line_plot(data[:,[0, state+1]],  x_units_in=['energy', 'a.u.'], y_units_in=['energy', 'au'], ax=myax, show_plot=False, save_plot=False)
        myax = myplot.line_plot(data[:,[0, nstates]],  x_units_in=['energy', 'a.u.'], y_units_in=['energy', 'au'], ax=myax, show_plot=True, save_plot=True)
        

if __name__ == "__main__":
    PESModel.from_inputfile('plot_pes.inp')
