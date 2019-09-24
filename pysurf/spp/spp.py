import configparser
import os
import numpy as np
import importlib


class SurfacePointProvider():
    """ The Surface Point Provider is the main interface providing the
        providing the information of a system at a specific point. It
        takes care where to take the information from according to the
        specified input file
    """
    def __init__(self, inputfile):
        """ The inputfile for the SPP has to provide the necessary
            information, how to produce the data at a specific point
            in the coordinate space.
        """
        self.config = configparser.ConfigParser()
        if os.path.isfile(inputfile):
            self.config.read(inputfile)

        """ If a model is used, import the model according to the user
            input and provide an instance in the variable self.user_inst
        """
        if self.config['MAIN']['mode'] == 'model':
            path = os.path.abspath(self.config['MODEL']['module'])
            spec = importlib.util.spec_from_file_location("", path)
            usr_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(usr_module)
            usr_class = getattr(usr_module, self.config['MODEL']['class'])
            self.usr_inst = usr_class()

    def get(self, coord):
        """ The get method is the method which should be called by
            external programs, which want to use the SPP. As an
            input it takes the coordinates and gives back the
            information at this specific position.
        """
        if self.config['MAIN']['mode'] == 'model':
            return self.usr_inst.get(coord)


if __name__ == "__main__":
    spp = SurfacePointProvider('./test.inp')
    print(spp.get(np.array([0.0, 0.0, 0.0])))
