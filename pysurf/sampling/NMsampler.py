from pysurf.colt Colt


class NMSampler(Colt):
    _questions = """
    number of points = 20 :: int

    # Range in which the normal modes should be sampled
    sampling range = [-5, 5] :: list

    # Decide whether sampling should be done along all normal modes
    selection of normal modes = all :: str :: [all, list]

    [selection of normal modes(all)]

    [selection of normal modes(list)]
    list of normal modes =  :: list

    moldenfile = :: str

    """

    def __init__(self, molecule, modes, is_massweighted=False):
        self.molecule = molecule
        self.modes = modes
        self.is_massweighted = is_massweighted
        self._check_modes()
        self.npoints = c


    def get_init(self):
        """Return all infos needed for the initial condition parser"""
        return {'molecule': self.molecule,
                'modes': self.modes}

    def get_condition(self):
        """Return a single created initial condition"""
        _, conds = get_initial_condition(self.molecule, self.modes)
        return conds

    @classmethod
    def from_config(cls, config):
        """ """
        if config['from'] == 'molden':
            return cls.from_molden(config['from']['moldenfile'])
        elif config['from'] == 'frequencies':
            return cls.from_freqs(config['from']['frequencies'])
        raise Exception("only (molden, frequencies) implemented")

    @classmethod
    def from_molden(cls, filename):
        molden = MoldenParser(filename, ['Info', 'Freqs', 'FrCoords', 'FrNormCoords'])
        # get molecule info
        atoms = [atom for atom, _, _, _ in molden['FrCoords']]
        atomids = np.array([atomname_to_id[atom] for atom in atoms])
        crd = np.array([[x, y, z] for _, x, y, z in molden['FrCoords']])
        masses = np.array([MASSES[idx]*U_TO_AMU for idx in atomids])
        # create molecule
        molecule = Molecule(atomids, crd, masses)
        #
        modes = [Mode(freq * CM_TO_HARTREE, np.array(molden['FrNormCoords'][imode]))
                 for imode, freq in enumerate(molden['Freqs'])]
        #
        modes = nm.create_mass_weighted_normal_modes(modes, molecule)
        #
        return cls(molecule, modes, True)

    def _check_modes(self):
        img = [mode.freq for mode in self.modes if mode.freq < 0.0]
        nimg_freq = len(img)
        if nimg_freq == 0:
            return

        def to_strg(number):
            return "%12.8f" % number

        print(f"Found {nimg_freq} imaginary frequencies:")
        print("[" + ", ".join(map(to_strg, img)) + "]")


       

