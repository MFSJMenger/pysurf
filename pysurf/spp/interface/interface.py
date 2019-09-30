from .qchem import QChem


class Interface():
    """ This interface class is just a in between class which
        checks which QM interface has to be called. 
        It should be extended when more QM programs are added.
    """
    def __init__(self, config, logger, refgeo, **kwargs):
        self.config = config
        self.logger = logger
        self.refgeo = refgeo
        self.atoms = refgeo['atoms']
        pass

    def get(self, coord):
        """ Get is the method that is called by the Surface Point
            Provider to start a QM calculation at a specific point
            and to get back the output.
        """
        if 'program' not in self.config.keys():
            self.logger.error('Program is not provided in '
                              + 'AB INITIO section!')
            exit()
        if self.config['program'] == 'qchem':
            if 'template' not in self.config.keys():
                self.logger.error('template not specified '
                                  + 'for QChem calculation!')
            qm = QChem(self.config, self.refgeo)
        else:
            self.logger.error('QM program not yet implemented')
            exit()
        result = qm.get(coord)
        if type(result) is dict:
            return result
        else:
            self.logger.error('Error in QM calculation: '+str(result))
