import os
import matplotlib.pyplot as plt
import numpy as np

from pysurf.sampling import Sampling
from pysurf.colt import Colt, FromCommandline
from pysurf.qctools.converter import Converter, time_converter
from pysurf.utils import exists_and_isfile
from pysurf.analysis import Plot
from pysurf.utils import SubfolderHandle
from pysurf.dynamics import DynDB

class NoFurtherQuestions(Colt):
    _questions = ""



class AnalysePopulation(Colt):

    _questions="""
    time_units = fs :: str :: [au, fs]
    save_data = yes :: str :: [yes, no]
    plot_population = yes :: str :: [yes, no]
    folder = prop :: str
    subfolder = traj :: str
    """

    _save_data = {
            'yes': "data_file = prop/population.dat :: file",
            'no' : " "
            }

    _plot_population = {
            'yes': "plot_inputfile = prop/plot_population.inp :: file", 
            'no' : ""
            }
    
    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases('save_data', {name: mode for name, mode in cls._save_data.items() })
        questions.generate_cases('plot_population', {name: mode for name, mode in cls._plot_population.items()})

    @classmethod
    def from_inputfile(cls, inputfile):
        if not(exists_and_isfile(inputfile)):
            config = cls.generate_input(inputfile, config=None)
        else:
            config = cls.generate_input(inputfile, config=inputfile)
        quests = cls.generate_questions(config=inputfile)
        config = quests.check_only(inputfile)
        

        if config['plot_population'].value == 'yes':
            plot_config = {}
        #    plot_config['x_label'] = 'time'
            plot_config['x_units'] = config['time_units']
            plot_config['x_label_unit'] = True
            plot_config['y_label_unit'] = False
            plot_config['y_label'] = 'population'
            plot_config = {'':plot_config}

            presets = ''
            presets += f"y_start = 0.0\n"
            presets += f"y_end = 1.0\n"
            presets += f"x_label = time\n"
            presets += f"y_label_unit = False\n"
            presets += f"y_label = population\n"
            presets += f"x_units = {config['time_units']}\n"
            presets += f"x_label_unit = True\n"
            presets += f"[save_plot(yes)]\nplot_file = prop/population.png\n"
            presets += "legend = {}\n"

            print(plot_config)
            plot_input = config['plot_population']['plot_inputfile']
            if exists_and_isfile(plot_input):
                plot_config = Plot.generate_input(plot_input, config=plot_input, presets=presets)
            else:
                plot_config = Plot.generate_input(plot_input, config=plot_config, presets=presets)
            print(plot_config)
        return cls(config, plot_config)


    def __init__(self, config, plot_config=None):
        self.config = config
        self.folder = config['folder']
        self.subfolder = config['subfolder']
        subfolderhandle = SubfolderHandle(self.folder, self.subfolder)
        propfiles = subfolderhandle.fileiter('prop.db')
        for idx, propfile in enumerate(propfiles):
            db = DynDB.load_database(propfile, read_only=True)
            dbtime = np.array(db['time']).flatten()
            if idx == 0:
                nsteps = len(db)
                nstates = db.info['dimensions']['nstates']
                data = np.zeros(shape=(nsteps, nstates + 1), dtype=float)
                data[:, 0] = dbtime
                counter = np.zeros(nsteps, dtype=int)
            if np.max(np.abs(data[:len(dbtime),0] - dbtime)) < 0.1:
                currstate = np.array(db['currstate']).flatten()
                for idx in range(len(dbtime)):
                    data[idx, int(currstate[idx]+1)] += 1
                    counter[idx] += 1
            else:
                print('times do not fit')
        
        data[:, 1:] = data[:, 1:]/np.repeat(counter, nstates).reshape(nsteps, nstates)

        converter = time_converter.get_converter('au', self.config['time_units'])
        data[:, 0] = converter(data[:, 0])

        myplt = Plot(plot_config)
        myax = myplt.line_plot(data[:,[0,1]], ('time', self.config['time_units']), y_units_in=None, ax=None, show_plot=False, save_plot=False, line_props={'label': 'state 0'})
        for state in range(1, nstates):
            save = False
            plot = False
            if state == nstates-1:
                save = True
                plot = True
            myax = myplt.line_plot(data[:,[0, state+1]], ('time', self.config['time_units']), y_units_in=None, ax=myax, show_plot=plot, save_plot=save, line_props={'label': f"state {state}"})

        




@FromCommandline("""
inputfile = prop/analyse_population.inp :: file
""")
def command_analyse_population(inputfile):
        analyse_population = AnalysePopulation.from_inputfile(inputfile)

if __name__ == "__main__":
    command_analyse_population()
