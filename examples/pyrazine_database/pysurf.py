import sys
import configparser
import logging
import os

inputfile = sys.argv[0]

logger = logging.getLogger('pysurf')
logger.setLevel(logging.INFO)
handler = logging.FileHandler(filename='pysurf.log', mode='w')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - '       
                               + '%(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('PYSURF started')

config = configparser.ConfigParser()

# read inputfile
if os.path.isfile(inputfile):
    config.read(inputfile)
    inputfile_full = os.path.abspath(inputfile)
    path = os.path.dirname(inputfile_full)
    config['MAIN']['path'] = path
else:
    logger.error('No valid inputfile found!')
    exit()

#set logging level
if 'logging' in config['MAIN'].keys():
    loglevel = config['MAIN']['logging']
    if loglevel == 'debug':
        logger.setLevel(logging.DEBUG)
        for hand in logger.handlers:
            hand.setLevel(logging.DEBUG)
        logger.info('logging level set to debug')
    else:
        logger.setLevel(logging.INFO)
        logger.info('logging level set to info as specification was'
                    + ' not known')

# get simulation
if 'simulation' in config['MAIN'].keys():
    if config['MAIN']['simulation'] == 'surface hopping':
        start_sh(config)
    else:
        self.logger.error('Simulation not known, stop program!')
        exit()

