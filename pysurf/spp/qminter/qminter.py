import os
from .qchem import QChem
import importlib


def get_qminter(config, logger, refgeo):
    if 'program' not in config.keys():
        logger.error('Program is not provided in '
                          + 'AB INITIO section!')
        exit()
    if config['program'] == 'qchem':
        if 'template' not in config.keys():
            logger.error('template not specified '
                              + 'for QChem calculation!')
        interface = QChem(config, refgeo)

    else:
        logger.error('QM program not yet implemented')
        exit()
    return interface
