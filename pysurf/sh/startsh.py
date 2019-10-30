import logging
import os
import numpy as np
import logging

from collections import namedtuple

from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.utils.constants import fs2au
from pysurf.utils.strutils import split_str
from pysurf.sh.landau_zener import landau_zener_surfacehopping
from pysurf.wigner import WignerSampling
from pysurf.wigner import InitialConditions

class StartSh():
    def __init__(self, config):
        self.logger = logging.getLogger('pysurf')

        # parse SURFACE HOPPING section
        if 'SURFACE HOPPING' in config.keys():
            try:
                self.restart = config.getboolean('SURFACE HOPPING', 'restart')
            except:
                self.logger.info('No information on restart, assume no restart!')
                self.restart = False
            


            # get whether specific trajectory is 

            # get number of steps
            try:
                self.nsteps = config.getint('SURFACE HOPPING', 'number of steps')
            except:
                self.logger.error('Could not get number of steps')
                exit()
            # get random seed for random number generator of hopping algorithm
            # only important for developpers and debugging
            try:
                self.seed = config.getint('SURFACE HOPPING', 'random seed')
            except:
                self.seed = None

            # get stepsize
            try:
                self.dt = config.getfloat('SURFACE HOPPING', 'stepsize')
                self.dt = self.dt * fs2au
            except:
                self.logger.info('Could not find stepsize. Assume 0.5 fs')
                self.dt = 0.5 * fs2au
            
            # get starting state
            try:
                self.initstate = config.getint('SURFACE HOPPING', 'initial state')
            except:
                if self.restart is False:
                    self.logger.error('No initial state is provided')
                    exit()

            # get number of trajectories
            try:
                self.ntrajs = config.getint('SURFACE HOPPING', 'number of trajectories')
            except:
                self.logger.info('Could not find number of trajectories')

            # get number of trajectory to run
            try:
                self.run_traj = config.get('SURFACE HOPPING', 'run trajectory')
                if self.run_traj.lower() == 'none' or self.run_traj.lower() == 'false':
                    self.run_traj = False
                else:
                    try: 
                        self.run_traj = int(self.run_traj)
                    except:
                        pass
            except:
                self.run_traj = None

            # get whether propagation is performed
            try:
                self.propagation = config.getboolean('SURFACE HOPPING', 'propagation')
            except:
                self.propagation = True

            # get inputfile for SPP
            try:
                self.sppinput = config.get('SURFACE HOPPING', 'inputfile for SPP')
                if os.path.isfile(self.sppinput):
                    self.sppinput = os.path.abspath(self.sppinput)
                else:
                    self.logger.error('Could not find inputfile for SPP')
            except:
                self.logger.info('SPP inputfile not given')
                exit()

            # get hopping algorithm
            try:
                self.hopalg = config.get('SURFACE HOPPING', 'hopping algorithm')
            except:
                self.logger.info('No hopping algorithm provided, assume landau zener')
                self.hopalg = 'landau zener'


            # get whether creation of initial conditions is needed
            try:
                self.preparation = config.getboolean('SURFACE HOPPING', 'preparation')
            except:
                self.preparation = True

        else:
            self.logger.error('No section SURFACE HOPPING in inputfile!')
            exit()

        # parse INITIAL CONDITIONS section
        self.init = False
        self.molden = False
        self.freqs = False
        self.icondb = False
        if ('INITIAL CONDITIONS' in config.keys()) and not(self.restart):

            # take initial conditions from file
            try:
                self.icondbfile = config.get('INITIAL CONDITIONS', 'from database file')
                if os.path.isfile(self.icondbfile):
                    self.icondbfile = os.path.abspath(self.icondbfile) 
                    self.preparation = False
                    self.icondb = True
                    self.logger.info('Using initial conditions from ' + str(self.icondbfile))
                else:
                    self.logger.error('Filepath of DB file with initial conditions is not correct')
                    exit()
                self.init = True
            except:
                pass

            try:
                self.sampling = config.get('INITIAL CONDITIONS', 'sampling')
            except:
                if self.preparation is True:
                    self.sampling = 'wigner'
                    self.logger.info('No sampling information found in INITIAL CONDITIONS! '
                                     + 'Assume Wigner sampling!')

            # molden input for initial conditions
            try:
                self.molden = config.get('INITIAL CONDITIONS', 'from molden file')
                if os.path.isfile(self.molden):
                    self.molden = os.path.abspath(self.molden)
                else:
                    self.logger.error('Molden file path not correct for INITIAL CONDITIONS')
                self.init = True
            except:
                pass

            #get frequencies for initial conditions
            try:
                vfloat = np.vectorize(float)
                self.freqs = vfloat(split_str(config.get('INITIAL CONDITIONS', 'frequencies')))
                self.init = True
            except:
                pass

        elif not('INITIAL CONDITIONS' in config.keys()) and not(self.restart):
            self.logger.error('No INITIAL CONDITIONS section in inputfile!')
        
        #decide whether to prepare initial conditions
        if self.preparation is True:
            initconds = self.create_init()
        elif self.icondb is True:
            initconds = self.read_init()

        # decide whether to run which trajectory
        if self.restart is False:
            self.start_sh(initconds)
        if self.restart is True:
            self.restart_sh()
    
    def create_init(self, filename='init_conds.db'):
        if self.init is False:
            self.logger.error('Not enough information provided for the initial sampling!')
            exit()
        if self.sampling == 'wigner':
            if self.molden is not False:
                sampling = WignerSampling.from_molden(self.molden)
                conditions = sampling.create_initial_conditions(filename, self.ntrajs)
            elif self.freqs is not False:
                sampling = WignerSampling.from_freqs(self.freqs)
                conditions = sampling.create_initial_conditions(filename, self.ntrajs, model=True)
            else:
                self.logger.error('No source for the sampling given, e.g. from molden file!')
        else:
            self.logger.error('Sampling method not known!')
            exit()
        return conditions

    def read_init(self):
        conditions = InitialConditions.from_db(self.icondbfile)
        return conditions

    def start_sh(self, initconds):
        print('Johannes start_sh: ', self.run_traj)
        if self.run_traj is None:
            rg = range(self.ntrajs)
        elif isinstance(self.run_traj, int):
            rg = range(self.run_traj, self.run_traj+1)
        else:
            rg = self.run_traj

        for i in rg:
            dirname = 'traj.' + str(i)
            if os.path.isdir(dirname):
                self.logger.info('Directory for trajectory ' + str(i) + ' already exists')
                continue
            os.mkdir(dirname)
            os.chdir(dirname)

            #setup sh logger
            shlogger = setup_shlogger()
            shlogger.info('starting trajectory ' + str(i))
            self.write_init_cond(initconds.get_condition(i))
            print('Johannes wrote initcond')
            print('Johannes self.propagation: ', self.propagation)
            if self.propagation is True:
                if self.hopalg == 'landau zener':
                    landau_zener_surfacehopping(initconds.get_condition(i),
                                                self.initstate,
                                                self.nsteps,
                                                self.seed,
                                                self.sppinput,
                                                self.dt)
                else:
                    self.logger.error('Hopping algorithm not implemented!')
                    exit()

            os.chdir('../')
            self.logger.info('Finished trajectory ' + str(i))
    
    def restart_sh(self):
        print('Johannes restart_sh: ', self.run_traj)
        if self.run_traj is None:
            rg = range(self.ntrajs)
        elif isinstance(self.run_traj, int):
            rg = range(self.run_traj, self.run_traj+1)
        else:
            rg = self.run_traj

        for i in rg:
            dirname = 'traj.' + str(i)
            if os.path.isdir(dirname):
                self.logger.info('Directory for trajectory ' + str(i) + ' already exists')
                os.chdir(dirname)

                #set up sh logger
                shlogger = setup_shlogger()

                if os.path.isfile('prop.db'):
                    shlogger.info('taking restart information from prop.db')
                    db = Database.load_db('prop.db')
                    nsdone = len(db['coord'])
                    if nsdone < self.nsteps:
                        init_cond = InitialCondition = namedtuple("InitialCondition", ['crd', 'veloc'])
                        init_cond.crd = db['coord'][-1]
                        init_cond.veloc = db['veloc'][-1]
                        self.initstate = db['curr_state'][-1][0]
                        shlogger.info('Restarting from step: ' + str(nsdone))
                        shlogger.info('Restarting from state: ' + str(self.initstate))
                        landau_zener_surfacehopping(init_cond,
                                                    self.initstate,
                                                    self.nsteps-nsdone,
                                                    self.seed,
                                                    self.sppinput,
                                                    self.dt)
                elif os.path.isfile('init_cond.db'):
                    shlogger.info('Taking restart information from init_cond.db')
                    self.icondbfile = 'init_cond.db'
 #                   print('Johannes in restart and looking for init_cond.db')
                    initcond = self.read_init_cond('init_cond.db')
                    landau_zener_surfacehopping(initcond,
                                        self.initstate,
                                        self.nsteps,
                                        self.seed,
                                        self.sppinput,
                                        self.dt)
                os.chdir('../')

    def write_init_cond(self, init_cond, filename='init_cond.db'):
        InitialConditions.save_condition(filename, init_cond)

    def read_init_cond(self, filename):
        return InitialConditions.from_db(filename).get_condition(0)

def setup_shlogger(level=logging.INFO):
    logger = logging.getLogger('sh')
    logger.setLevel(level)

    for hand in logger.handlers:
        hand.stream.close()
        logger.removeHandler(hand)
    handler = logging.FileHandler(filename='sh.log', mode='a')
    handler.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - '       
                                  + '%(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger
