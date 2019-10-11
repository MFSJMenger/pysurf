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
        if 'coord' in request.keys():
            coord = request['coord']
        else:
            for prop in self.properties:
                request[prop] = None
            return request

        if 'interpolation' in self.config.keys():
            try:
                interpol = self.config.getboolean('interpolation')
            except ValueError:
                interpol = self.config['interpolation']
        else:
            interpol = False

        if interpol is True:
            interpol = 'linear'
            self.logger.info('set default interpolation to linear'
                             + 'interpolation')

        self.logger.debug('interpolation method is ' + str(interpol))
        if interpol:
            # if db empty, start qm calculation
            if len(self.db['coord']) == 0:
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

            if 'max distance interpolation' in self.config.keys():
                thr = float(self.config['max distance interpolation'])
            else:
                thr = 0.25
                self.logger.info('Using default max distance for'
                                 + 'interpolation of 0.25')

            if min_dist['dist'] < thr:
                res = self.interpolator.get(coord)
                success = True
                for key in res.keys():
                    if res[key] is None:
                        success = False
                if success is False:
                    self.logger.debug('Interpolation was not'
                                      + 'successfull! QM'
                                      + 'calculation is started!')
                    self.get_qm(request)
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
                self.energy = LinearNDInterpolator(coords, energies)
                self.gradient = LinearNDInterpolator(coords, gradients)
                self.logger.info('successfully set up interpolator')
            except:
                self.energy = lambda x: None
                self.gradient = lambda x: None
        else:
            self.logger.error('ERROR: Interpolation scheme not yet'
                              + 'implemented!')
            exit()

    def get(self, coord):
        self.logger.debug('using interpolator to get energy')
        en = self.energy(coord)[0]
        grad = self.gradient(coord).reshape((self.nstates,
                                             self.natoms, 3))

        return {'energy': en, 'gradient': grad, 'coord': coord}
