"""
PySurf Module:
    Validation and Training of Interpolators

Provide infrastructure for the training of interpolators
and test them against a validation set
"""
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
from pysurf.colt import from_commandline
from pysurf.analysis import Plot3D






class AnalyseFitPes2Nm(Colt):
    _questions = """
        spp = spp.inp :: existing_file
        mode =  :: int
        mode2 = :: int
        energy_units = eV :: str :: [eV, au, cm-1, nm]
        reference_energy = 0 :: float
        plot_input = plot_fit_pes_2nm.inp :: file
        start = -2 :: float
        end = 2 :: float
        start2 = -2 :: float
        end2 = 2 :: float
        npoints = 15 :: int
        npoints2 = 15 :: int
        save_data = yes :: str :: [yes, no]
        plot_pes = yes :: str :: [yes, no]
        moldenfile = molden.in :: existing_file
        states = :: ilist, optional
    """

    _save_data = {
            'yes': "data_file = fit_pes_2nm.dat :: file",
            'no' : " "
            }

    _plot_pes = {
            'yes': "plot_inputfile = plot_fit_pes_2nm.inp :: file", 
            'no' : ""
            }


    @classmethod
    def _extend_questions(cls, questions):
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
            plot_config['z_label_unit'] = True
            plot_config['z_label'] = 'energy'
            plot_config = {'':plot_config}

            presets = ''
            presets += f"z_label = energy\n"
            presets += f"x_label_unit = True\n"
            presets += f"[save_plot(yes)]\nplot_file = plot_fit_pes_nm.png\n"
            presets += "legend = {}\n"

            plot_input = config['plot_pes']['plot_inputfile']
            if exists_and_isfile(plot_input):
                plot_config = Plot3D.generate_input(plot_input, config=plot_input, presets=presets)
            else:
                plot_config = Plot3D.generate_input(plot_input, config=plot_config, presets=presets)
            print(plot_config)
        return cls(config, plot_config)
        

    def __init__(self, config, plot_config):
        #
        self.logger = get_logger('plotnm.log', 'plotnm', [])
        #
        #Get Surface Point Provider
        config_spp = self._get_spp_config(config['spp'])
        natoms, self.nstates, properties = self._get_db_info(config_spp['use_db']['database'])
        atomids = [1 for _ in range(natoms)]
        self.spp = SurfacePointProvider.from_config(config_spp, properties, self.nstates, natoms, atomids=atomids,
                                        logger=self.logger)
        #
        self.interpolator = self.spp.interpolator
        self.interpolator.train()
        #
        qx, qy, crds = self.generate_crds(config['moldenfile'], mode=(config['mode'], config['mode2']), start=(config['start'], config['start2']), end=(config['end'], config['end2']), npoints=(config['npoints'], config['npoints2']))
        energy = self._compute(crds)
        data = []
        conv = energy_converter.get_converter('au', config['energy_units'])
        for qxi, qyi, en in zip(np.ravel(qx), np.ravel(qy), energy):
            en = conv(en) - conv(config['reference_energy'])
            data += [[qxi, qyi, *en]]
        data = np.array(data)
        nstates = data.shape[1] - 2
        energy = data[:,2:]
        energy = energy.T.reshape((nstates, *qx.shape))
        
        # Take only states according to user input
        if config['states'] is None:
            statelist = [i for i in range(nstates)]
        else:
            statelist = config['states']

        if config['save_data'] == 'yes':
            np.savetxt(config['save_data']['data_file'], data)

        if config['plot_pes'] == 'yes':
            myplt = Plot3D(plot_config)

            save = False
            plot = False
            for state in statelist:
                if state == nstates-1:
                    save = True
                    plot = True
                if state == statelist[0]:
                    myax = myplt.surface_plot((qx, qy, energy[state]), 
                                       x_units_in=('length', 'au'), 
                                       y_units_in=('length', 'au'), 
                                       z_units_in=('energy', config['energy_units']), 
                                       ax=None, 
                                       show_plot=plot, 
                                       save_plot=save)
                else:
                    myax = myplt.surface_plot((qx, qy, energy[state]), 
                                       x_units_in=('length', 'au'), 
                                       y_units_in=('length', 'au'), 
                                       z_units_in=('energy', config['energy_units']), 
                                       ax=myax, 
                                       show_plot=plot, 
                                       save_plot=save)



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

        crds = np.empty((npoints[0]*npoints[1], *crd.shape))
        qx, qy = np.mgrid[start[0]: end[0]: npoints[0]*1j, start[1]:end[1]:npoints[1]*1j]
        i=0
        for qxi, qyi in zip(np.ravel(qx), np.ravel(qy)):
            crds[i] = crd + qxi*modes[mode[0]].displacements + qyi*modes[mode[1]].displacements
            i += 1
        return qx, qy, crds


    def _compute(self, crds, states=None):
        ndata = len(crds)
        energy = []
        for i, crd in enumerate(crds):
            if states is None:
                result = self.spp.request(crd, ['energy'])
            else:
                result = self.spp.request(crd, ['energy'])
            #
            energy += [np.copy(result['energy'])]
        return np.array(energy)

@from_commandline("""
inputfile = analyse_fit_pes_2nm.inp :: file
""")
def command_analyse_pes(inputfile):
        analyse_pes = AnalyseFitPes2Nm.from_inputfile(inputfile)


if __name__ == '__main__':
    command_analyse_pes()
