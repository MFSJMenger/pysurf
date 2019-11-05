import numpy as np
from scipy.interpolate import LinearNDInterpolator
from copy import deepcopy

from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.utils.strutils import split_str
from ..qminter.qminter import get_qminter


class DBInter():
    """This class handels all the interaction with the database and
        the surface point provider, 
        saves the data and does the interpolation
    """
    def __init__(self, dbpath, config, logger, refgeo):
        self.dbpath = dbpath
        self.config = config
        self.logger = logger
        self.refgeo = refgeo
        self.natoms = len(self.refgeo['atoms'])
        self.nstates = int(self.config['number of states'])

        try:
            self.inter_grad = self.config.getboolean('interpolate gradient')
            if self.inter_grad is None:
                self.inter_grad = False
                self.config['interpolate gradient'] = 'False'
        except:
            self.inter_grad = False
            self.config['interpolate gradient'] = 'False'


        if 'properties' in config.keys():
            self.properties = split_str(config['properties'])
        else:
            self.properties = ['energy', 'gradient']
        self.db = self.create_db(self.dbpath)
        if 'interpolation' in self.config.keys():
            if len(self.db['coord']) > 0:
                self.interpolator = Interpolation(self.db, self.config,
                                                  self.logger,
                                                  self.refgeo)

        if 'interpolation' in self.config.keys():
            try:
                self.interpol = self.config.getboolean('interpolation')
            except ValueError:
                self.interpol = self.config['interpolation']
        else:
            self.interpol = False
        if self.interpol is True:
            self.interpol = 'linear'
            self.logger.info('set default interpolation to linear'
                             + 'interpolation')
        self.logger.debug('interpolation method is ' + str(self.interpol))

        if self.interpol:
            if 'max distance interpolation' in self.config.keys():
                self.thr = float(self.config['max distance interpolation'])
            else:
                self.thr = 0.25
                self.logger.info('Using default max distance for'
                                 + ' interpolation of 0.25')


    def create_db(self, filename='db.dat'):
        variables = {}
        variables['coord'] = DBVariable(np.double, ('frame', 'natoms',
                                                    'three'))
        if 'energy' in self.properties:
            variables['energy'] = DBVariable(np.double, ('frame',
                                                         'nstates'))
        if 'gradient' in self.properties:
            variables['gradient'] = DBVariable(np.double, ('frame',
                                                           'nstates',
                                                           'natoms',
                                                           'three'))


        dct = {'dimensions': {
                              'frame': 'unlimited',
                              'natoms': self.natoms,
                              'nstates': self.nstates,
                              'three': 3,
                              'one': 1},
               'variables': variables
               }
        db = Database(filename, dct)
        return db

    def get(self, request):
        self.logger.debug('DBInter.get is called')
        if 'coord' in request.keys():
            coord = request['coord']
        else:
            for prop in self.properties:
                request[prop] = None
            return request


        if self.interpol:
            self.logger.debug('Number of DB entries: ' 
                               + str(len(self.db['coord'])))
            # if db empty, start qm calculation
            if len(self.db['coord']) == 0:
                self.logger.debug('DB is empty, start qm calc')
                res = self.get_qm(request)
                return res

            # get closest geometry
            min_dist = self.get_min_dist(coord)
            self.logger.debug('smallest distance of new coord is: '
                              + str(min_dist['dist']) + ' of entry: '
                              + str(min_dist['entry']))

            if min_dist['dist'] < self.thr:
                # If point is already in DB
                if min_dist['dist'] < 0.001:
                    self.logger.info('take point from DB: entry '
                                     + str(min_dist['entry']))
                    for key in self.properties:
                        request[key] = self.db[key][min_dist['entry']]
                    if self.inter_grad is True:
                        request['gradient'] = self.interpolator.calc_grad(coord)
                    return request
                
                res = self.interpolator.get(coord)
                success = True
                for key in res.keys():
                    if res[key] is None:
                        success = False
                if success is False:
                    self.logger.debug('Interpolation was not'
                                      + ' successfull! QM'
                                      + ' calculation is started!')
                    res = self.get_qm(request)

            else:
                res = self.get_qm(request)

        else:
            res = self.get_qm(request)

        return res

    def get_qm(self, request):
        qminter = get_qminter(self.config, self.logger, self.refgeo)
        res = qminter.get(request)
        increase = False
        if 'coord' in res.keys():
            self.db.append('coord', res['coord'])
            increase = True
        if self.inter_grad is False:
            if 'gradient' in res.keys():
                self.db.append('gradient', res['gradient'])
        if 'energy' in res.keys():
            self.db.append('energy', res['energy'])
        if increase: self.db.increase
        
        self.interpolator = Interpolation(self.db, self.config,
                                                  self.logger,
                                                  self.refgeo)
        if self.inter_grad is True:
            res['gradient'] = self.interpolator.calc_grad(request['coord'])

        return res

    def get_min_dist(self, coord, method='max atom displacement'):
        if method == 'max atom displacement':
            for i in range(len(self.db['coord'])):
                max_dist = np.amax(np.absolute(coord
                                   - self.db['coord'][i]))
                if i == 0:
                    min_dist = max_dist
                    res = {'entry': i, 'coord': self.db['coord'][i],
                           'dist': min_dist}
                elif max_dist < min_dist:
                    min_dist = max_dist
                    res = {'entry': i, 'coord': self.db['coord'][i],
                           'dist': min_dist}
        return res


class Interpolation():
    def __init__(self, db, config, logger, refgeo):
        self.db = db
        self.config = config
        self.logger = logger
        self.natoms = len(refgeo['atoms'])
        self.nstates = int(self.config['number of states'])

        # The interpolate gradient key has been added in DBInter
        self.inter_grad = self.config.getboolean('interpolate gradient')
        
        if self.config['interpolation'] == 'linear':
            noc = len(self.db['coord'])
            coords = np.empty((noc, self.natoms*3), dtype=float)
            energies = np.empty((noc, self.nstates), dtype=float)
            gradients = np.empty((noc, self.nstates*self.natoms*3),
                                 dtype=float)
            for i in range(len(self.db['coord'])):
                coords[i, :] = self.db['coord'][i].flatten()
                energies[i] = self.db['energy'][i]
                gradients[i, :] = self.db['gradient'][i].flatten()
            try:
                self.logger.info('try to set up interpolator')
                self.energy = LinearNDInterpolator(coords, energies)
                self.gradient = LinearNDInterpolator(coords, gradients)
                self.logger.info('successfully set up interpolator')
            except:
                #self.energy = LinearNDInterpolator(coords, energies)
                self.energy = lambda x: None
                self.gradient = lambda x: None

        if self.config['interpolation'] == 'shepard':
            self.logger.info('Using Shepard interpolation scheme')
            noc = len(self.db['coord'])
            coords = np.empty((noc, self.natoms*3), dtype=float)
            energies = np.empty((noc, self.nstates), dtype=float)
            if self.inter_grad is False:
                gradients = np.empty((noc, self.nstates*self.natoms*3),
                                 dtype=float)
            for i in range(len(self.db['coord'])):
                coords[i, :] = self.db['coord'][i].flatten()
                energies[i] = self.db['energy'][i]
                if self.inter_grad is False:
                    gradients[i, :] = self.db['gradient'][i].flatten()
            self.logger.info('setting up Shepard interpolator')
            self.energy = ShepardInterpolator(coords, energies)

            #Check whether gradient is calculated or interpolated

            if self.inter_grad is False:
                self.gradient = ShepardInterpolator(coords, gradients)
        else:
            self.logger.error('ERROR: Interpolation scheme not yet'
                              + 'implemented!')
            exit()

    def get(self, coord):
        self.logger.debug('using interpolator to get energy')
        en = self.energy(coord.flatten())
        if self.inter_grad is False:
            grad = self.gradient(coord.flatten())
        else:
            grad = self.calc_grad(coord)
        grad.resize((self.nstates,*coord.shape))
#        grad.resize((self.nstates,len(coord),3))
#         try:
#             en = self.energy(coord)[0]
#             print('Johannes: ', en)
#         except:
#             en = None
#         try:
#             grad = self.gradient(coord).reshape((self.nstates,
#                                              self.natoms, 3))
#         except:
#             grad = None
        return {'energy': en, 'gradient': grad, 'coord': coord}

    def calc_grad(self, coord, dq=0.01):
        grad = np.empty((self.nstates, len(coord.flatten())), dtype=float)

        for i in range(len(coord.flatten())):
            coord1 = deepcopy(coord).flatten()
            coord1[i] += dq
            en1 = self.energy(coord1)
            coord2 = deepcopy(coord).flatten()
            coord2 -= dq
            en2 = self.energy(coord2)
            grad[:,i] = (en1 - en2)/2.0/dq

        grad.resize((self.nstates, *coord.shape))
        return grad


class ShepardInterpolator():
    def __init__(self, coords, values):
        self.coords = coords
        self.values = values
    
    def __call__(self, coord):
        weights = self.get_weights(coord)
        res = 0.0
        for i, value in enumerate(self.values):
            res += weights[i]*value
        res = res/sum(weights)
        return res

    def get_weights(self,coord):
        weights = np.zeros(len(self.coords), dtype=float)
        for i in range(len(self.coords)):
            weights[i] = 1./np.linalg.norm((coord-self.coords[i]))**2
        return weights
