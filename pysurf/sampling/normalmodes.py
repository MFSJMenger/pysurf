from copy import deepcopy
#
from collections import namedtuple
#
import numpy as np
#
from ..constants import U_TO_AMU
#

Mode = namedtuple("Mode", ["freq", "displacements"])


class NormalModes(object):

    @staticmethod
    def is_normal_mode_format(modes, natoms):
        """Check if the normal mode is in the given format!
            change the docstring to something useful!
        """

        thresh = 0.05
        nmodes = len(modes)
        #
        matrix = np.array([[mode.displacements[iatom][ixyz] for mode in modes]
                          for iatom in range(natoms)
                          for ixyz in range(3)])
        #
        result = np.dot(matrix.T, matrix)  # returns modes*modes matrix
        # compute the trace
        trace = np.trace(result)
        # set diagonal elements to 0.0
        np.fill_diagonal(result, 0.0)
        #
        nmodes_m1 = float(nmodes-1)
        # check that trace == nmodes - 1
        if trace > nmodes_m1:
            if abs((trace/nmodes) - 1.0) < thresh:
                # compute max/min ofdiagonal elements
                max_val = abs(np.max(result))
                min_val = abs(np.min(result))
                # check that there are no significant elements!
                if max_val < thresh and min_val < thresh:
                    return True
        return False

    @staticmethod
    def compute_norm_mode(mode, molecule):
        """Compute the norm of a mode"""
        norm = 0.0
        for iatom, displacement in enumerate(mode.displacements):
            for xyz in displacement:
                norm += xyz**2 * molecule.masses[iatom]/U_TO_AMU
        return np.sqrt(norm)

    @staticmethod
    def scale_mode(imode, modes, expression):
        """example expression:

        lambda imode, iatom, ixyz, disp: disp/(norm/np.sqrt(molecule.masses[iatom]))

        """
        mode = deepcopy(modes[imode])
        #
        for iatom, displacement in enumerate(mode.displacements):
            mode.displacements[iatom] = [expression(imode, iatom, ixyz, disp)
                                         for ixyz, disp in enumerate(displacement)]
        #
        return mode

    @classmethod
    def create_mass_weighted_normal_modes(cls, modes, molecule):
        """decide which format the normal modes are and convert to the mass-weighted once"""
        ANG_TO_BOHR = 1./0.529177211
        # compute norms
        norms = [cls.compute_norm_mode(mode, molecule) for mode in modes]
        # define options
        options = {
                'gaussian-type': lambda imode, iatom, ixyz, disp: (disp/(norms[imode]
                                                                   / np.sqrt(molecule.masses[iatom]
                                                                   / U_TO_AMU))),
                'molpro-type': lambda imode, iatom, ixyz, disp: (disp * np.sqrt(
                                                                 molecule.masses[iatom]
                                                                 / U_TO_AMU)),
                'columbus-type': lambda imode, iatom, ixyz, disp: (disp * np.sqrt(
                                                                   (molecule.masses[iatom]
                                                                    / U_TO_AMU)
                                                                   / ANG_TO_BOHR)),
                'mass-weighted': lambda imode, iatom, ixyz, disp: disp,
        }
        #
        possible_formats = {}
        #
        for typ, expression in options.items():
            current_modes = [cls.scale_mode(imode, modes, expression)
                             for imode in range(len(modes))]
            if cls.is_normal_mode_format(current_modes, molecule.natoms):
                possible_formats[typ] = current_modes
        #
        if len(possible_formats) < 1:
            print("possible_formats = ", possible_formats)
            raise Exception('Could not specify format possible formats = "%s"'
                            % ", ".join(possible_formats.keys()))
        #
        if len(possible_formats) > 1:
            print("possible_formats = ", ", ".join(possible_formats.keys()))
            print("select first one...")
        #
        for _, mode in possible_formats.items():
            return mode

    @classmethod
    def create_dimensionless_normal_modes(cls, modes, molecules, mass_weighted=False):
        if mass_weighted is False:
            # get mass weighted normal modes
            modes = cls.create_mass_weighted_normal_modes(modes, molecules)
        modes = [Mode(mode.freq, mode.displacements/np.sqrt(mode.freq)) for mode in modes]
        return modes
