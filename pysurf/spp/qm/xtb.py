import os
import subprocess as sp
#
from . import AbinitioBase
from ...fileparser import  read_geom
from .xtbhelp import XTBReader
from ...system import Molecule


def run_interface(script, name):
    print('############## Run section: ##############\nInterface %s call:' % name)
    print(script)
    error=sp.call(script, shell=True)
    if error==0:
        print('Finished!')
    else:
        print('*** Something went wrong! ***')
        print("Error count = %d" % error)
    return


def write_content_to_file(fileName, content, options="w"):
    """ write content to file [fileName]

        fileName = str, Name of the file to be read
        content  = str, content written to the file
    """
    with open(fileName, options) as f:
        f.write(content)


class XTBInterface(AbinitioBase):

    name = "XTB"

    _user_input = """
    executable = xtb
    name = xtb calculation
    refgeom = :: existing_file
    filename = xtb_input.xyz
    """

    implemented = ['energy', 'gradient']

    def __init__(self, executable, name, atomids, filename):
        self.natoms = len(atomids)
        self.molecule = Molecule(atomids, None)
        self.name = name
        self.executable = executable
#        self.natoms, self.atomnames, self.coords = read_geom(refgeom)
        self.ifile = filename

    @classmethod
    def from_config(cls, config, atomids, nstates, nghost):
        assert nghost == 0
        if nstates > 1:
            raise Exception("Only Groundstate calculations possible")
        return cls(config['executable'], config['name'], atomids, config['filename'])
    
    def get(self, request):
        self.molecule.crd = request.crd
        self._write_general_inputs()
        (en, dip), grad, _ = self._run_gradient()
        if 'gradient' in request:
            request['gradient'][0] = grad.reshape((self.natoms, 3))
        request.set('energy', en)
        return request

    def _write_general_inputs(self):
        self.molecule.write_xyz(self.ifile)

    def _run_energy(self):
        command = "%s %s --sp --copy > output_energy" % (self.executable, self.ifile)
        run_interface(command, self.name)
        return (self._read_xtb_energy_and_dipole("output_energy"), None, None)

    def _run_gradient(self):
        command = "%s %s --grad --copy > output_gradient " % (self.executable,  self.ifile)
        run_interface(command, self.name)
        output = (self._read_xtb_energy_and_dipole("output_gradient"), self._read_xtb_gradient("gradient"), None)
        self.clean_up()
        return output

    def _run_frequency(self):
        (_, _), grad, _ = self._run_gradient()
        command = "%s %s --hess --copy > output_freq" % (self.executable, self.ifile)
        run_interface(command, self.name)
        return (self._read_xtb_energy_and_dipole("output_freq"), 
                grad, 
                self._read_xtb_hessian("hessian"))


    def clean_up(self):
        os.rename("gradient", "grad_old")
        return

    def _read_xtb_energy_and_dipole(self, filename):
        xtb = XTBReader(filename, ["Dipole", "Energy"])
        return xtb["Energy"], xtb["Dipole"]

    def _read_xtb_gradient(self, filename):
        xtb = XTBReader(filename, ["Gradient"], {"NAtoms": self.natoms})
        return xtb["Gradient"]

    def _read_xtb_hessian(self, filename):
        xtb = XTBReader(filename, ["Hessian"], {"NAtoms": self.natoms})
        return xtb["Hessian"]

    @classmethod
    def get_coordinates(cls, atomid, coords):
        return "% 4s    % 14.10f   % 14.10f   % 14.10f \n"  % (atomid, coords[0], coords[1], coords[2])
