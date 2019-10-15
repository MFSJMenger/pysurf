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

    # for test purposes also the model input is possible
    elif config['program'] == 'model':
        path = config['path']
        try: 
            # get absolute path for model 
            path = os.path.join(path, config['module']) 
            path = os.path.realpath(path) 
            logger.info('The model is given as: ' + path) 
            spec = importlib.util.spec_from_file_location("", path) 
            usr_module = importlib.util.module_from_spec(spec) 
            spec.loader.exec_module(usr_module) 
        except (FileNotFoundError, ModuleNotFoundError, ImportError): 
            logger.error('The user module could not be ' 
                              + 'loaded: ' + path) 
            exit() 
        try: 
            usr_class = getattr(usr_module, config['class']) 
        except AttributeError: 
            logger.error('The user class could not be found: ' 
                              + config['MODEL']['class']) 
            exit() 
        interface = usr_class()
    else:
        logger.error('QM program not yet implemented')
        exit()
    return interface
