import numpy as np
from scipy.interpolate import LinearNDInterpolator

from pysurf.database.database import Database
from pysurf.database.dbtools import DBVariable
from pysurf.utils.strutils import split_str
from ..qminter.qminter import get_qminter


class DBInter():
    """This class handels all the interaction with the database and
        the surface point provider, i.e. saves the stuff and does the
        interpolation
    """
    def __init__(self, dbpath, config, logger, refgeo):
        self.dbpath = dbpath
        self.config = config
        self.logger = logger
        self.refgeo = refgeo
        self.natoms = len(self.refgeo['atoms'])
        self.nstates = int(self.config['number of states'])
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
                if len(self.db['coord']) > 0:
                    self.interpolator = Interpolation(self.db,
                                                      self.config,
                                                      self.logger,
                                                      self.refgeo)
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
                    self.interpolator = Interpolation(self.db,
                                                      self.config,
                                                      self.logger,
                                                      self.refgeo)
            else:
                res = self.get_qm(request)
                self.interpolator = Interpolation(self.db, self.config,
                                                  self.logger,
                                                  self.refgeo)
        else:
            res = self.get_qm(request)

        print('Johannes in dbinter.get() res= ', res)
        return res

    def get_qm(self, request):
        qminter = get_qminter(self.config, self.logger, self.refgeo)
        res = qminter.get(request)
        increase = False
        if 'coord' in res.keys():
            self.db.append('coord', res['coord'])
            increase = True
        if 'gradient' in res.keys():
            self.db.append('gradient', res['gradient'])
        if 'energy' in res.keys():
            self.db.append('energy', res['energy'])
        if increase: self.db.increase
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
                print('Johannes energies:', energies)
                print('Johannes grads:', gradients)
                print('Johannes coords:', coords)
                self.energy = LinearNDInterpolator(coords, energies)
                print('Johannes energies set up')
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
            gradients = np.empty((noc, self.nstates*self.natoms*3),
                                 dtype=float)
            for i in range(len(self.db['coord'])):
                coords[i, :] = self.db['coord'][i].flatten()
                energies[i] = self.db['energy'][i]
                gradients[i, :] = self.db['gradient'][i].flatten()
            self.logger.info('setting up Shepard interpolator')
            self.energy = ShepardInterpolator(coords, energies)
            self.gradient = ShepardInterpolator(coords, gradients)
        else:
            self.logger.error('ERROR: Interpolation scheme not yet'
                              + 'implemented!')
            exit()

    def get(self, coord):
        self.logger.debug('using interpolator to get energy')
        en = self.energy(coord.flatten())
        grad = self.gradient(coord.flatten())
        print('Johannes: ', en)
        print('Johannes: ', grad)
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

class ShepardInterpolator():
    def __init__(self, coords, values):
        self.coords = coords
        self.values = values
    
    def __call__(self,coord):
        print('Johannes called Shepard interpolator')
        weights = self.get_weights(coord)
        res = 0
        for i in range(len(self.values)):
            res += weights[i]*self.values[i]
        res = res/sum(weights)
        return res

    def get_weights(self,coord):
        weights = np.zeros(len(self.coords), dtype=float)
        for i in range(len(self.coords)):
            weights[i] = np.linalg.norm((coord-self.coords[i]))
        return weights