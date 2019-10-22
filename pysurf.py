# -*- coding: utf-8 -*-

"""Main module."""

import sys
import os
import logging
import configparser

from pysurf.sh.startsh import StartSh

logger = logging.getLogger('pysurf')
logger.setLevel(logging.INFO)
handler = logging.FileHandler(filename='pysurf.log', mode='w')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - '       
                              + '%(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('PYSURF')

config = configparser.ConfigParser()

if len(sys.argv) > 1:
    inputfile = sys.argv[1]
    if os.path.isfile(inputfile):
        config.read(inputfile)
        inputfile_full = os.path.abspath(inputfile)
        path = os.path.dirname(inputfile_full)
        config['MAIN']['path'] = path
    else:
        logger.error('Inputfile ' + inputfile + ' for pysurf not found!')
        exit()

    if 'logging' in config['MAIN'].keys():
        loglevel = config['MAIN']['logging']
        if loglevel == 'debug':
            logger.setLevel(logging.DEBUG)
            for hand in logger.handlers:
                hand.setLevel(logging.DEBUG)
            logger.info('logging level set to debug')
        else:
            logger.info('logging level set to info as specification'
                        + ' was not known')

    if 'simulation' in config['MAIN'].keys():
        if config['MAIN']['simulation'] == 'surface hopping':
            simulation = StartSh(config)

else:
    logger.error('Inputfile for pysurf not given!')
    exit()

