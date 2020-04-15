import os
import shutil
import numpy as np
import yt
import scipy.constants as scc
import glob

"""
This file is part of the suite of scripts to use LibEnsemble on top of WarpX
simulations. It reads the output plotfiles of a simulation and returns beam
quantity relevant for the LibEnsemble optimizations ('f') as well as other
beam quantities for convenience.
"""

yt.funcs.mylog.setLevel(50)


def _beam_properties(filepath):
    """
    Reads plotfile filepath and compute and return beam parameters

    Parameters
    ----------
    filepath : path to plotfile to read
    """

    # Read beam quantities from plotfile
    ds = yt.load(filepath)
    ad = ds.all_data()
    w = ad['beam', 'particle_weight'].v
    x = ad['beam', 'particle_position_x'].v
    ux = ad['beam', 'particle_momentum_x'].v/scc.m_e/scc.c
    uy = ad['beam', 'particle_momentum_y'].v/scc.m_e/scc.c
    uz = ad['beam', 'particle_momentum_z'].v/scc.m_e/scc.c

    # Compute beam parameters
    # Defined like that, the beam charge is > 0.
    charge = np.sum(w) * scc.e
    gamma = np.sqrt(1. + ux**2 + uy**2 + uz**2)
    energy_MeV = .511 * (gamma - 1.)
    energy_avg = np.mean(energy_MeV)
    energy_std = np.std(energy_MeV) / energy_avg
    emittance = np.sqrt(np.std(x)**2 * np.std(ux)**2 - np.mean(x*ux)**2)
    return charge, energy_avg, energy_std, emittance


def read_sim_output(workdir):
    """
    Return optimizing quantity 'f' and other parameters for convenience.

    Parameters
    ----------
    workdir : Path to directory where the simulation ran.
    """
    # Get beam properties at the beginning of the run
    datafile = 'diags/plotfiles/plt00000/'
    filepath = os.path.join(workdir, datafile)
    charge_i, _, _, emittance_i = _beam_properties(filepath)

    # Get beam properties at the end of the run
    file_list = glob.glob('diags/plotfiles/plt?????')
    file_list.sort()
    datafile = file_list[-1]
    filepath = os.path.join(workdir, datafile)
    charge_f, energy_avg, energy_std, emittance_f = _beam_properties(filepath)

    # delete simulation results, just to have smaller data
    shutil.rmtree('diags')

    # Build a quantity to minimize (f) that encompasses
    # emittance AND charge loss 1% charge loss has the
    # same impact as doubling the initial emittance.
    # we minimize f!
    f = emittance_f + emittance_i*(1.-charge_f/charge_i)*100
    warpx_out = np.array([f, energy_std, energy_avg, charge_f, emittance_f])

    return warpx_out
