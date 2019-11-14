import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import Rbf
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
            if len(self.db['crd']) > 0:
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
            self.interpol = 'rbf'
            self.logger.info('set default interpolation to rbf'
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
        variables['crd'] = DBVariable(np.double, ('frame', 'natoms',
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
        if 'crd' in request.keys():
            crd = request['crd']
        else:
            for prop in self.properties:
                request[prop] = None
            return request


        if self.interpol:
            self.logger.debug('Number of DB entries: ' 
                               + str(len(self.db['crd'])))
            # if db empty, start qm calculation
            if len(self.db['crd']) == 0:
                self.logger.debug('DB is empty, start qm calc')
                res = self.get_qm(request)
                return res

            # get closest geometry
            min_dist = self.get_min_dist(crd)
            self.logger.debug('smallest distance of new crd is: '
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
                        request['gradient'] = self.interpolator.calc_grad(crd)
                    return request
                
                res = self.interpolator.get(crd)
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
            
        print('Johannes res:', res)
        return res

    def get_qm(self, request):
        qminter = get_qminter(self.config, self.logger, self.refgeo)
        res = qminter.get(request)
        increase = False
        if 'crd' in res.keys():
            self.db.append('crd', res['crd'])
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
            res['gradient'] = self.interpolator.calc_grad(request['crd'])

        return res

    def get_min_dist(self, crd, method='max atom displacement'):
        if method == 'max atom displacement':
            for i in range(len(self.db['crd'])):
                max_dist = np.amax(np.absolute(crd
                                   - self.db['crd'][i]))
                if i == 0:
                    min_dist = max_dist
                    res = {'entry': i, 'crd': self.db['crd'][i],
                           'dist': min_dist}
                elif max_dist < min_dist:
                    min_dist = max_dist
                    res = {'entry': i, 'crd': self.db['crd'][i],
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
        
        if self.config['interpolation'] == 'rbf':
            noc = len(self.db['crd'])
            crds = np.empty((noc, self.natoms*3), dtype=float)
            energies = np.empty((noc, self.nstates), dtype=float)
            if self.inter_grad is False:
                gradients = np.empty((noc, self.nstates*self.natoms*3),
                                 dtype=float)
            for i in range(len(self.db['crd'])):
                crds[i, :] = self.db['crd'][i].flatten()
                energies[i] = self.db['energy'][i]
                if self.inter_grad is False:
                    gradients[i, :] = self.db['gradient'][i].flatten()
            try:
                self.logger.info('try to set up interpolator')
                self.energy = MyRBF(crds, energies)
                if self.inter_grad is False:
                    self.gradient = MyRBF(crds, gradients)
                self.logger.info('successfully set up interpolator')
            except:
                self.energy = lambda x: None
                if self.inter_grad is False:
                    self.gradient = lambda x: None
                else:
                    self.gradient = lambda x: np.zeros(slef.nstates, self.atoms, 3, dtype=float)
        elif self.config['interpolation'] == 'shepard':
            self.logger.info('Using Shepard interpolation scheme')
            noc = len(self.db['crd'])
            crds = np.empty((noc, self.natoms*3), dtype=float)
            energies = np.empty((noc, self.nstates), dtype=float)
            if self.inter_grad is False:
                gradients = np.empty((noc, self.nstates*self.natoms*3),
                                 dtype=float)
            for i in range(len(self.db['crd'])):
                crds[i, :] = self.db['crd'][i].flatten()
                energies[i] = self.db['energy'][i]
                if self.inter_grad is False:
                    gradients[i, :] = self.db['gradient'][i].flatten()
            self.logger.info('setting up Shepard interpolator')
            self.energy = ShepardInterpolator(crds, energies)

            #Check whether gradient is calculated or interpolated

            if self.inter_grad is False:
                self.gradient = ShepardInterpolator(crds, gradients)
        else:
            self.logger.error('ERROR: Interpolation scheme not yet'
                              + 'implemented!')
            exit()

    def get(self, crd):
        self.logger.debug('using interpolator to get energy')
        en = self.energy(crd.flatten())
        if self.inter_grad is False:
            print('Johannes: interpolate gradient')
            grad = self.gradient(crd.flatten())
        else:
            print('Johannes: calculated gradient')
            grad = self.calc_grad(crd)
        grad.resize((self.nstates,*crd.shape))
#        grad.resize((self.nstates,len(crd),3))
#         try:
#             en = self.energy(crd)[0]
#         except:
#             en = None
#         try:
#             grad = self.gradient(crd).reshape((self.nstates,
#                                              self.natoms, 3))
#         except:
#             grad = None
        return {'energy': en, 'gradient': grad, 'crd': crd}

    def calc_grad(self, crd, dq=0.01):
        grad = np.zeros((self.nstates, len(crd.flatten())), dtype=float)

        if self.energy(crd) is None:
            pass
        else:
            for i in range(len(crd.flatten())):
                crd1 = deepcopy(crd).flatten()
                crd1[i] += dq
                en1 = self.energy(crd1)
                crd2 = deepcopy(crd).flatten()
                crd2[i] -= dq
                en2 = self.energy(crd2)
                grad[:,i] = (en1 - en2)/2.0/dq

        grad.resize((self.nstates, *crd.shape))
        return grad



class ShepardInterpolator():
    def __init__(self, crds, values):
        self.crds = crds
        self.values = values
    
    def __call__(self, crd):
        weights = self._get_weights(crd)
        res = 0.0
        for i, value in enumerate(self.values):
            res += weights[i]*value
        res = res/sum(weights)
        return res

    def _get_weights(self,crd):
        weights = np.zeros(len(self.crds), dtype=float)
        for i in range(len(self.crds)):
            weights[i] = 1./np.linalg.norm((crd-self.crds[i]))**2
        return weights


class MyRBF():
    def __init__(self, crds, values):
        print('Setting up rbf')
        if len(crds.shape) == 1:
            crds = deepcopy(crds)
            crds.resize(len(crds),1)
        if len(values.shape) > 1:
            self.vdim = len(values[0])
            self.rbfi = []
            for i in range(self.vdim):
                self.rbfi += [Rbf(*crds.transpose(), values[:,i])]
        else:
            self.vdim = 1
            self.rbfi = Rbf(*crds.transpose(), values)

    def __call__(self, crd):
#        crd = deepcopy(crdin)
        res = []
        if self.vdim != 1:
            for i in range(self.vdim):
                res += [self.rbfi[i](*crd.flatten())]
        else:
            res = self.rbfi(crd)
        return np.array(res).flatten()
