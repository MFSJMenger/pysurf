"""
PySurf Module:
    Validation and Training of Interpolators

Provide infrastructure for the training of interpolators
and test them against a validation set
"""
import numpy as np

from pysurf.database import PySurfDB
from pysurf.spp import SurfacePointProvider
from pysurf.logger import get_logger
from pysurf.colt import Colt
from pysurf.sampling.normalmodes import Mode
from pysurf.system import Molecule
from pysurf.sampling.normalmodes import NormalModes as nm
from pysurf.molden import MoldenParser
from pysurf.constants import U_TO_AMU, CM_TO_HARTREE
from pysurf.system.atominfo import ATOMNAME_TO_ID, MASSES
from pysurf.utils import exists_and_isfile
from pysurf.qctools.converter import energy_converter
from pysurf.colt import FromCommandline

from pysurf.analysis import Plot

from scipy.optimize import minimize




class AnalyseFitPesNm(Colt):

    _questions = """
        spp = spp.inp :: existing_file
        savefile = :: str
        mode =  :: int
        energy_units = eV :: str :: [eV, au, cm-1, nm]
        reference_energy = 0 :: float
        plot_input = plot_fit_pes_nm.inp :: file
        start = -2 :: float
        end = 2 :: float
        npoints = 100 :: int
        save_data = yes :: str :: [yes, no]
        plot_pes = yes :: str :: [yes, no]
        moldenfile = molden.in :: existing_file
    """

    _save_data = {
            'yes': "data_file = fit_pes_nm.dat :: file",
            'no' : " "
            }

    _plot_pes = {
            'yes': "plot_inputfile = plot_fit_pes_nm.inp :: file", 
            'no' : ""
            }


    @classmethod
    def _generate_subquestions(cls, questions):
        questions.generate_cases('save_data', {name: mode for name, mode in cls._save_data.items() })
        questions.generate_cases('plot_pes', {name: mode for name, mode in cls._plot_pes.items()})


    @classmethod
    def from_inputfile(cls, inputfile):
        if not(exists_and_isfile(inputfile)):
            config = cls.generate_input(inputfile, config=None)
        else:
            config = cls.generate_input(inputfile, config=inputfile)
        quests = cls.generate_questions(config=inputfile)
        config = quests.check_only(inputfile)
        
        plot_config={}
        if config['plot_pes'].value == 'yes':
            plot_config = {}
            plot_config['x_label'] = 'mode'
            plot_config['x_label_unit'] = False
            plot_config['y_label_unit'] = True
            plot_config['y_label'] = 'energy'
            plot_config = {'':plot_config}

            presets = ''
            presets += f"y_label = energy\n"
            presets += f"x_label_unit = True\n"
            presets += f"[save_plot(yes)]\nplot_file = plot_fit_pes_nm.png\n"
            presets += "legend = {}\n"

            plot_input = config['plot_pes']['plot_inputfile']
            if exists_and_isfile(plot_input):
                plot_config = Plot.generate_input(plot_input, config=plot_input, presets=presets)
            else:
                plot_config = Plot.generate_input(plot_input, config=plot_config, presets=presets)
        return cls(config, plot_config)
        

    def __init__(self, config, plot_config):
        #
        self.logger = get_logger('plotnm.log', 'plotnm', [])
        #
        #Get Surface Point Provider
        config_spp = self._get_spp_config(config['spp'])
        natoms, self.nstates, properties = self._get_db_info(config_spp['use_db']['database'])
        atomids = [1 for _ in range(natoms)]
        self.spp = SurfacePointProvider(None, properties, self.nstates, natoms, atomids,
                                        logger=self.logger, config=config_spp)
        #
        self.interpolator = self.spp.interpolator
        self.savefile = config['savefile']
        self.interpolator.train(self.savefile)
        #
        qs, crds = self.generate_crds(config['moldenfile'], mode=config['mode'], start=config['start'], end=config['end'], npoints=config['npoints'])
        energy = self._compute(crds)
        data = []
        conv = energy_converter.get_converter('au', config['energy_units'])
        for q, en in zip(qs, energy):
            en = conv(en) - conv(config['reference_energy'])
            data += [[q, *en]]
        data = np.array(data)
        
        if config['save_data'] == 'yes':
            np.savetxt(config['save_data']['data_file'], data)
        
        if config['plot_pes'] == 'yes':
            nstates = data.shape[1]-1
            myplt = Plot(plot_config)
            myax = myplt.line_plot(data[:,[0,1]], x_units_in=('length', 'au'), y_units_in=('energy', 'au'), ax=None, show_plot=False, save_plot=False)
            for state in range(1, nstates):
                save = False
                plot = False
                if state == nstates-1:
                    save = True
                    plot = True
                myax = myplt.line_plot(data[:,[0, state+1]], x_units_in=('length', 'au'), y_units_in=('energy', 'au'), ax=myax, show_plot=plot, save_plot=save)


    def _get_spp_config(self, filename):
        questions = SurfacePointProvider.generate_questions(presets="""
                use_db=yes :: yes
                [use_db(yes)]
                fit_only = yes :: yes
                write_only = no :: no
                """)
        return questions.ask(config=filename, raise_read_error=False)

    def _get_db_info(self, database):
        db = PySurfDB.load_database(database, read_only=True)
        rep = db.dbrep
        natoms = rep.dimensions.get('natoms', None)
        if natoms is None:
            natoms = rep.dimensions['nmodes']
        nstates = rep.dimensions['nstates']
        return natoms, nstates, db.saved_properties

    def generate_crds(self, moldenfile, mode, start=-5, end=5, npoints=50):
        molden = MoldenParser(moldenfile, ['Info', 'Freqs', 'FrCoords', 'FrNormCoords'])
        # get molecule info
        atoms = [atom for atom, _, _, _ in molden['FrCoords']]
        atomids = np.array([ATOMNAME_TO_ID[atom] for atom in atoms])
        crd = np.array([[x, y, z] for _, x, y, z in molden['FrCoords']])
        masses = np.array([MASSES[idx]*U_TO_AMU for idx in atomids])
        # create molecule
        molecule = Molecule(atomids, crd, masses)
        #
        print('nmsampler, molecule', molecule)
        modes = [Mode(freq * CM_TO_HARTREE, np.array(molden['FrNormCoords'][imode]))
                 for imode, freq in enumerate(molden['Freqs'])]
        #
        modes = nm.create_mass_weighted_normal_modes(modes, molecule)

        crds = np.empty((npoints, *crd.shape))
        qs = np.linspace(start, end, npoints)
        for i, q in enumerate(qs):
            crds[i] = crd + q*modes[mode].displacements
        return qs, crds



    def validate(self, filename, properties):
        db = PySurfDB.load_database(filename, read_only=True)
        self._compute(db, properties)

    def save_pes(self, filename, database):
        db = PySurfDB.load_database(database, read_only=True)
        results, _ = self._compute(db, ['energy'])

        def str_join(values):
            return ' '.join(str(val) for val in values)

        with open(filename, 'w') as f:
            f.write("\n".join(f"{i} {str_join(fitted)} {str_join(exact)}" 
                              for i, (fitted, exact) in enumerate(results['energy'])))

    def _compute(self, crds):
        ndata = len(crds)
        energy = []
        for i, crd in enumerate(crds):
            result = self.spp.request(crd, ['energy'])
            #
            energy += [np.copy(result['energy'])]
        return energy

@FromCommandline("""
inputfile = analyse_fit_pes_nm.inp :: file
""")
def command_analyse_pes(inputfile):
        analyse_pes = AnalyseFitPesNm.from_inputfile(inputfile)


if __name__ == '__main__':
    command_analyse_pes()
