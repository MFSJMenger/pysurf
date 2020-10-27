import numpy as np

from ..spp import SurfacePointProvider
from ..database import PySurfDB
from ..sampling import Sampler

from . import engine

@engine.register_action
def spp_analysis(sppinp: "file") -> "spp":
    """ comment 
        
        Parameters
        ----------
            -sppinp: file
                blabla

        Returns
        -------
            spp
                blabla4
    """
    #
    config = _get_spp_config_analysis(sppinp)
    #
    natoms, nstates, properties = _get_db_info(config['use_db']['database'])
    atomids = [1 for _ in range(natoms)]
    spp = SurfacePointProvider.from_config(config, properties, nstates, natoms, atomids,
                                    logger=None)
    #
    return spp

@engine.register_action
def spp_calc(sppinp: "file", atomids: "ilist", nstates: "int", properties: "list"=['energy']) -> "spp":
    #
    presets=f"""
            use_db=yes :: yes
            [use_db(yes)]
            write_only = yes :: yes
            properties = {properties} :: {properties}
            """
    config = SurfacePointProvider.generate_input(sppinp, config=sppinp, presets=presets)
    natoms = len(atomids)
    spp = SurfacePointProvider.from_config(config, properties, nstates, natoms, nghost_states=0, atomids=atomids,
                                    logger=None)
    #
    return spp

@engine.register_action
def get_energies(spp: "spp", crds: "crds") -> "array2D":
    energies = []
    for crd in crds:
        request = spp.request(crd, ['energy'])
        energies += [request['energy']] 
    return np.array(energies)

@engine.register_action
def sampler(samplerinp: "file") -> "sampler":
    sampler = Sampler.from_inputfile(samplerinp)
    return sampler

@engine.register_action
def crds_from_sampler(sampler: "sampler", npoints: "int") -> "crds":
    crds = []
    for i in range(npoints):
        cond = next(sampler)
        crds += [cond.crd]
    return np.array(crds)

@engine.register_action
def sp_calc(spp: "spp", crd: "crd", properties: "list"=['energy']) -> "request":
    return spp.request(crd, properties)

def _get_spp_config_analysis(filename):
    questions = SurfacePointProvider.generate_questions(presets="""
            use_db=yes :: yes
            [use_db(yes)]
            write_only = no :: no
            [use_db(yes)::interpolator]
            fit_only = yes :: yes
            """)
    return questions.ask(config=filename, raise_read_error=False)

def _get_db_info(database):
    db = PySurfDB.load_database(database, read_only=True)
    rep = db.dbrep
    natoms = rep.dimensions.get('natoms', None)
    if natoms is None:
        natoms = rep.dimensions['nmodes']
    nstates = rep.dimensions['nstates']
    return natoms, nstates, db.saved_properties

    
