import os
import matplotlib.pyplot as plt
import numpy as np

from pysurf.sampling import Sampling
from pysurf.colt import Colt, from_commandline 
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
    folder = ./ :: str
    subfolder = traj :: str
    nsteps = :: int, optional
    """

    _save_data = {
            'yes': "data_file = population.dat :: file",
            'no' : " "
            }

    _plot_population = {
            'yes': "plot_inputfile = plot_population.inp :: file", 
            'no' : ""
            }
    
    @classmethod
    def _extend_questions(cls, questions):
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
        
        plot_config={}
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
            presets += f"[save_plot(yes)]\nplot_file = population.png\n"
            presets += "legend = {}\n"

            plot_input = config['plot_population']['plot_inputfile']
            if exists_and_isfile(plot_input):
                plot_config = Plot.generate_input(plot_input, config=plot_input, presets=presets)
            else:
                plot_config = Plot.generate_input(plot_input, config=plot_config, presets=presets)
        return cls(config, plot_config)


    def __init__(self, config, plot_config=None):
        self.config = config
        self.folder = os.path.abspath(config['folder'])
        self.subfolder = config['subfolder']
        subfolderhandle = SubfolderHandle(self.folder, self.subfolder)
        propfiles = subfolderhandle.fileiter('prop.db')
        first = True
        for idx, propfile in enumerate(propfiles):
            db = DynDB.load_database(propfile, read_only=True)
            print(propfile)
            dbtime = np.array(db['time']).flatten()
            if len(dbtime) == 0:
                continue
            if first is True:
                first = False
                if config['nsteps'] is not None:
                    nsteps = config['nsteps']
                else:
                    nsteps = len(db)
                nstates = db.info['dimensions']['nstates']
                data = np.zeros(shape=(nsteps, nstates + 1), dtype=float)
                data[:min(len(db), nsteps), 0] = dbtime
                counter = np.zeros(nsteps, dtype=int)
                if nsteps > len(db):
                    times_missing = True
                    time_steps = len(db)
                else:
                    times_missing = False
                    time_steps = len(db)
            if np.max(np.abs(data[:min(time_steps, len(dbtime)),0] - dbtime[:min(time_steps, len(dbtime))])) < 0.1:
                if len(dbtime) > time_steps and times_missing is True:
                    data[time_steps:len(dbtime), 0] = dbtime[time_steps:min(nsteps, len(dbtime))]
                    time_steps = len(dbtime)
                currstate = np.array(db['currstate']).flatten()
                for idx in range(min(len(dbtime), nsteps)):
                    data[idx, int(currstate[idx]+1)] += 1
                    counter[idx] += 1
            else:
                print('times do not fit')
        
        data[:, 1:] = data[:, 1:]/np.repeat(counter, nstates).reshape(nsteps, nstates)

        converter = time_converter.get_converter('au', self.config['time_units'])
        data[:, 0] = converter(data[:, 0])
        
        if config['save_data'].value == 'yes':
            np.savetxt(config['save_data']['data_file'], data)

        if config['plot_population'].value == 'yes':
            myplt = Plot(plot_config)
            myax = myplt.line_plot(data[:,[0,1]], ('time', self.config['time_units']), y_units_in=None, ax=None, show_plot=False, save_plot=False, line_props={'label': 'state 0'})
            for state in range(1, nstates):
                save = False
                plot = False
                if state == nstates-1:
                    save = True
                    plot = True
                myax = myplt.line_plot(data[:,[0, state+1]], ('time', self.config['time_units']), y_units_in=None, ax=myax, show_plot=plot, save_plot=save, line_props={'label': f"state {state}"})

        




@from_commandline("""
inputfile = analyse_population.inp :: file
""")
def command_analyse_population(inputfile):
        analyse_population = AnalysePopulation.from_inputfile(inputfile)

if __name__ == "__main__":
    command_analyse_population()
