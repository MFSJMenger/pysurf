Running Nonadiabatic Trajectories with pysurf_fssh and SAOOVQE
==============================================================

This guide describes how to set up and run nonadiabatic molecular dynamics trajectories
using the `pysurf_fssh` and `SAOOVQE` repositories.

---

1. Environment Setup
---------------------

Clone the repositories:

    git clone https://github.com/qc2nl/SAOOVQE.git
    git clone https://github.com/gatox/pysurf_fssh.git

We recommend using a Conda environment to manage dependencies. First, create and activate a new environment:

    conda create -n pysurf_saoovqe python=3.9
    conda activate pysurf_saoovqe

Install required dependencies for SAOOVQE:

    pip install numpy scipy numba openfermion cirq
    conda install -c psi4 psi4

Then install the SAOOVQE package (note that this must be done **before** setting up PySurf):

    cd SAOOVQE
    pip install -e .

Install dependencies for PySurf manually, as it does not support pip installation:

    pip install pycolt qctools netcdf4 numpy matplotlib

---

2. Environment Variables
-------------------------

Add the following to your `~/.bash_profile` or `~/.zshrc` (depending on your shell):

    export PYTHONPATH=/full/path/to/pysurf_fssh:$PYTHONPATH
    export SAOOVQEDIR=/full/path/to/SAOOVQE

Make sure to replace `/full/path/to/...` with the actual paths on your machine.

---

3. Recommended Branches
------------------------

Use the following branches for compatibility:

    - pysurf_fssh: old_branch
    - SAOOVQE: noise_old_saoovqe

Switch branches as follows:

    cd pysurf_fssh
    git checkout old_branch
    cd ../SAOOVQE
    git checkout noise_old_saoovqe

---

4. Preparing Initial Conditions
-------------------------------

Create a folder for storing trajectory data:

    mkdir test_example_trajectories
    cd test_example_trajectories

Ensure the directory contains the Molden file (e.g., `CH2NH_molden`) for vibrational normal modes.

Run the sampling script:

    python /path/to/pysurf_fssh/bin/sampling.py

Example terminal inputs:

    n_conditions_max [100]: 10
    method (Wigner/NMSampler): Wigner
    sampling_db: <press Enter>
    [Wigner::from(molden)]
    moldenfile: CH2NH_molden

---

5. Setting Up the Trajectories
-------------------------------

    python /path/to/pysurf_fssh/bin/setup_propagation.py

The script will interactively request input via the terminal:

    n_traj [-1]: <press Enter>
    sampling_db: sampling.db
    spp: spp.inp
    initial state: 1
    prop: prop.inp
    copy_db: <press Enter>
    logging: <press Enter>
    mode: ab-initio
    software: INTSAOOVQE
    chg [0]: <press Enter>
    mult [1]: <press Enter>
    basis: cc-pvdz
    nelec_active: 4
    frozen_indices: 6
    active_indices: [6, 9]
    virtual_indices: [9, 43]
    do_oo_process: True
    noise: True
    mean: 0.0
    variance: 1e-8
    noise_energy_after_sa_vqe: True
    noise_before_orb_opt_phase: False
    noise_final_state_resolution: True
    noise_vqe_cost_function_energy: False
    noise_rdms_gradient: True
    noise_tdms_nacs: True
    properties: <press Enter>
    write_only: yes
    database: db.dat
    t: 0.0
    dt: 5
    mdsteps: 5
    substeps: false
    instate: 1
    nstates: 2
    states: [0, 1]
    ncoeff: [0.0, 1.0]
    prob: tully
    rescale_vel: nacs
    rev_vel_no_hop: True
    coupling: nacs
    method: Surface_Hopping
    decoherence: EDC
    n_substeps: False
    res_nacs: True

---

6. Notes
---------

- The Molden file must contain vibrational normal modes and can be generated using Q-Chem or OpenMolcas.
- Input files created (`spp.inp`, `prop.inp`, etc.) are editable after initial generation.

License: Apache License 2.0
