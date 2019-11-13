from string import Template
import numpy as np
import os

from pysurf.utils.osutils import exists_and_isfile
from pysurf.utils.strutils import split_str
from pysurf.utils.constants import bohr2angstrom


class QChem():
    def __init__(self, config, refgeo):
        """ class to write the input file, start the calculation and
            read the output for QChem.
            The class has to be initialized and then get is used 
            to start the calculation.
        """
        self.config = config
        self.refgeo = refgeo
        self.nstates = int(config['number of states'])
        self.ref_en = float(config['reference energy'])
        self.natoms = len(refgeo['atoms'])
        if 'properties' in config.keys():
            self.properties = split_str(config['properties'])
        else:
            self.properties = ['energy', 'gradient']
        try:
            self.method = config['method']
        except KeyError:
            self.method = 'dft'

    def get(self, request):
        """ get is the method that should be used by the interface.
            It takes as arguments the coordinates at a specific
            position and starts the QChem calculation there.
            Finally it reads the output and provides the data.
        """
        if 'crd' in request.keys():
            crd = request['crd']
            self.write_input(self.refgeo['atoms'], crd)
            outfile = self.start_calc()
            ret = self.read_output(outfile)
            if 'energy' in self.properties:
                request['energy'] = ret['energy']
            if 'gradient' in self.properties:
                request['gradient'] = ret['gradient']
            return request
        else:
            for prop in self.properties:
                request[prop] = None
            return request

    def write_input(self, atoms, crd, filename='qchem.in'):
        """ Prepares a input file from the template file.
            Here only the geometry is added.
        """
        self.config['template'] = os.path.join(self.config['path'],
                                               self.config['template'])
        with open(self.config['template']) as infile:
            temp = Template(infile.read())

        crdstring = ''
        for i in range(len(atoms)):
            crdstring += '{0} {1:12.5f} {2:12.5f} {3:12.5f}\n'\
                .format(atoms[i], crd[i, 0]*bohr2angstrom,
                        crd[i, 1]*bohr2angstrom,
                        crd[i, 2]*bohr2angstrom)
        with open(filename, 'w') as outfile:
            outfile.write(temp.safe_substitute(
                          geometry=crdstring.strip('\n')))

    def start_calc(self, filename='qchem.in'):
        """ Starts the QChem calculation by calling an external bash
            script so that everybody can adopt the bash script to 
            their cluster.
        """
#        path_inter = os.path.join(self.config['run command'])
        filename_input = filename
        if filename_input.split('.')[-1] == 'inp':
            filename_output = filename_input.strip('inp') + 'out'
        elif filename_input.split('.')[-1] == 'in':
            filename_output = filename_input.strip('in') + 'out'
        else:
            filename_output = filename_input + '.out'

        # If simulate is given in the inputfile, QChem calculation
        # will not be started
        simulate = False
        if 'simulate' in self.config.keys():
            if self.config.getboolean('simulate'):
                simulate = True
        if simulate is not True:
            # start QChem calculation
            cmd = self.config['run command']+' '+filename_input+' '+filename_output
            os.system(cmd)
        return filename_output

    def read_output(self, filename):
        """ Method to read the output file of the QChem calculation.
            At the moment energies and gradients are read from SCF 
            and TDDFT calculations.
        """
        if not(exists_and_isfile(filename)):
            return 'Outputfile of QChem calculation not found!'
            exit()
        with open(filename) as infile:
            res = infile.readlines()
        ret = {}
        if self.method == 'dft':
            if 'energy' in self.properties:
                ret['energy'] = self.read_energies_dft(res)
            if 'gradient' in self.properties:
                ret['gradient'] = self.read_gradients_dft(res)

        if self.method == 'adc':
            if 'energy' in self.properties:
                ret['energy'] = self.read_energies_adc(res)
            if 'gradient' in self.properties:
                ret['gradient'] = self.read_gradients_adc(res)

        return ret

    def read_energies_dft(self, res):
        en = np.zeros((self.nstates), dtype=float)
        en_check = np.empty_like(en, dtype=bool)
        for line in res:
            # Get ground state energy
            if 'energy in the final basis set' in line:
                en[0] = float(line.split('=')[-1]) - self.ref_en
                en_check[0] = True
            # Get excited state energy
            if 'Total energy for state' in line:
                state = int(line.split()[4].strip(':'))
                # sometimes QChem prints too many states
                if state <= self.nstates:
                    en[state] = float(line.split()[5]) - self.ref_en
                    en_check[state] = True
        # Check that energies of all states were found
        if np.all(en_check):
            return en
        else:
            return 'Not all energies were found in QChem output!'

    def read_gradients_dft(self, res):
        grad = np.zeros((self.nstates, self.natoms, 3), dtype=float)
        grad_check = np.zeros(self.nstates, dtype=bool)
        natoms = self.natoms
        vfloat = np.vectorize(float)
        for i in range(len(res)):
            if 'Gradient of' in res[i]:
                # Get the state of the gradient
                if 'Gradient of SCF Energy' in res[i]:
                    state = 0
                if 'Gradient of the state energy (including CIS Excitation '\
                        'Energy)' in res[i]:
                    state = int(res[i-1].split()[1])

                # Read in the gradient
                for j in range(int(natoms/6)):
                    for k in range(3):
                        grad[state, j*6:j*6+6, k] = vfloat(
                            res[i+j*4+2+k].split()[1:7])
                if natoms % 6 > 0:
                    j = int(natoms/6)
                    rest = natoms % 6
                    for k in range(3):
                        grad[state, j*6:j*6+rest, k] = vfloat(
                            res[i+j*4+2+k].split()[1:])
                grad_check[state] = True
        # Check that the gradients of all states were found
        if grad_check.all():
            return grad
        else:
            return 'Not all gradients were found in QChem output!'

    def read_energies_adc(self, res):
        en = np.zeros((self.nstates), dtype=float)
        en_check = np.empty_like(en, dtype=bool)
        for i in range(len(res)):
            # Get ground state energy
            if 'MP(2) Summary' in res[i]:
                en[0] = float(res[i+3].split()[-2]) - self.ref_en
                en_check[0] = True
            # Get excited state energy
            if 'Excited state' in res[i]:
                state = int(res[i].split()[2])
                en[state] = float(res[i+4].split()[2]) - self.ref_en
                en_check[state] = True
        # Check that energies of all states were found
        if np.all(en_check):
            return en
        else:
            return 'Not all energies were found in QChem output!'


    def read_gradients_adc(self, res):
        grad = np.zeros((self.nstates, self.natoms, 3), dtype=float)
        grad_check = np.zeros(self.nstates, dtype=bool)
        natoms = self.natoms
        vfloat = np.vectorize(float)
        for i in range(len(res)):
            # state for the gradient has to be found in a complicated
            # way. Here it is taken from the user input, which is shown
            # in the output of the calculation.
            if 'adc_state_to_opt' in res[i].lower():
                state = int(split_str(res[i].split()[-1])[-1])

            # Read actual gradient
            if 'Full Analytical Gradient' in res[i]:
                # Get the state of the gradient
                if 'Gradient of MP2 Energy' in res[i]:
                    state = 0

                # Read in the gradient
                for j in range(int(natoms/5)):
                    for k in range(3):
                        grad[state, j*5:j*5+5, k] = vfloat(
                            res[i+j*4+2+k].split()[1:7])
                if natoms % 5 > 0:
                    j = int(natoms/5)
                    rest = natoms % 5
                    for k in range(3):
                        grad[state, j*5:j*5+rest, k] = vfloat(
                            res[i+j*4+2+k].split()[1:])
                grad_check[state] = True
        # Check that the gradients of all states were found
        if grad_check.all():
            return grad
        else:
            return 'Not all gradients were found in QChem output!'

