import os
import matplotlib.pyplot as plt
import numpy as np

from pysurf.sampling import Sampling
from pysurf.colt import Colt, FromCommandline
from pysurf.qctools.converter import Converter, energy_converter
from pysurf.utils import exists_and_isfile
from pysurf.analysis import Plot

class NoFurtherQuestions(Colt):
    _questions = ""


class Broadening(Colt):
    _questions = """ 
        method = Lorentzian :: str :: [Gaussian, Lorentzian]
        width = :: float
        energy_start = :: float
        energy_end = :: float
        n_points = 500 :: int
    """

    def __init__(self, config):
        self.config = config

    def broadening(self, data):
        res = np.zeros((self.config['n_points'], 2))
        res[:,0] = np.linspace(self.config['energy_start'], self.config['energy_end'], self.config['n_points'])
        for  (e0, fosc) in data:
            if self.config['method'] == 'Lorentzian':
                func = self.l_vec(fosc, e0, self.config['width']) 
            if self.config['method'] == 'Gaussian':
                func = self.g_vec(fosc, e0, self.config['width'])
            res[:, 1] += func(res[:, 0])
        res[:, 1] = res[:, 1] / len(data)
        return res
    
    @staticmethod
    def l_vec(fosc, e0, w):
        return np.vectorize(lambda e:fosc * 1 / (1 + ((e - e0)/w*2)**2))

    @staticmethod
    def g_vec(fosc, e0, w):
        return np.vecotrize(lambda e: fosc * np.exp(-np.ln(2)*((e - e0)/w*2)**2))


class AnalyseSpectrum(Colt):
    specfile = os.path.join('spectrum', 'spectrum.db') 
    
    _questions="""
    input_db = spectrum/spectrum.db :: existing_file
    energy_units = au :: str :: [au, eV]
    broadening = yes :: str :: [yes, no]
    save_data = yes :: str :: [yes, no]
    plot_spectrum = yes :: str :: [yes, no]
    """

    _save_data = {
            'yes': "data_file = spectrum/spectrum.dat :: file",
            'no' : " "
            }

    _broadening={
            'yes': Broadening,
            'no' : NoFurtherQuestions
            }

    _plot_spectrum = {
            'yes': "plot_inputfile = spectrum/plot_spectrum.inp :: file", 
            'no' : ""
            }
    
    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases('save_data', {name: mode for name, mode in cls._save_data.items() })
        questions.generate_cases('broadening', {name: mode.questions for name, mode in cls._broadening.items()})
        questions.generate_cases('plot_spectrum', {name: mode for name, mode in cls._plot_spectrum.items()})


    @classmethod
    def from_inputfile(cls, inputfile):
        if not(exists_and_isfile(inputfile)):
            config = cls.generate_input(inputfile, config=None)
        else:
            config = cls.generate_input(inputfile, config=inputfile)
        quests = cls.generate_questions(config=inputfile)
        config = quests.check_only(inputfile)
        

        if config['plot_spectrum'].value == 'yes':
            plot_config = {}
            if config['broadening'].value == 'yes':
                plot_config['x_start'] = config['broadening']['energy_start']
                plot_config['x_end'] = config['broadening']['energy_end']
            plot_config['x_label'] = 'energy'
            plot_config['y_label_unit'] = True
            plot_config['y_label'] = 'intensity'
            plot_config['x_units'] = config['energy_units']
            plot_config['x_label_unit'] = True
            plot_config = {'': plot_config}

            presets = ''
            if config['broadening'].value == 'yes':
                presets += f"x_start = {config['broadening']['energy_start']}\n"
                presets += f"x_end = {config['broadening']['energy_end']}\n"
            presets += f"x_label = energy\n"
            presets += f"y_label_unit = True\n"
            presets += f"y_label = intensity\n"
            presets += f"y_units = a.u.\n"
            presets += f"x_units = {config['energy_units']}\n"
            presets += f"x_label_unit = True\n"
            presets += f"[save_plot(yes)]\nplot_file = spectrum/spectrum.png\n"

            plot_input = config['plot_spectrum']['plot_inputfile']
            if exists_and_isfile(plot_input):
                plot_config = Plot.generate_input(plot_input, config=plot_input, presets=presets)
            else:
                plot_config = Plot.generate_input(plot_input, config=plot_config, presets=presets)

        return cls(config, plot_config)


    def __init__(self, config, plot_config=None):
        self.config = config
        sampling = Sampling.from_db(self.specfile)

        nstates = sampling.info['dimensions']['nstates']
        npoints = sampling.nconditions
        data = []
        for  energy, fosc in  zip(np.copy(sampling._db['energy']), np.copy(sampling._db['fosc'])):
            for idx, en in enumerate(energy[1:]):
                data += [[en - energy[0], fosc[idx+1]]]
        data = np.array(data)

        converter = energy_converter.get_converter('au', self.config['energy_units'])
        data[:, 0] = converter(data[:, 0])
        if config['broadening'].value == 'yes':
            cont_data = Broadening(config['broadening']).broadening(data)
        else:
            cont_data = np.sort(data)

        Plot(plot_config).line_plot(cont_data, ('energy', self.config['energy_units']), y_units_in=None, ax=None, line_props={'linestyle':'solid'})
        




@FromCommandline("""
inputfile = spectrum/analyse_spectrum.inp :: file
""")
def command_analyse_spectrum(inputfile):
        analyse_spectrum = AnalyseSpectrum.from_inputfile(inputfile)

if __name__ == "__main__":
    command_analyse_spectrum()
