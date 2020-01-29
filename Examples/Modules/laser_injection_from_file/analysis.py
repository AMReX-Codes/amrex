#!/usr/bin/env python3

# Copyright 2020 Andrew Myers, Axel Huebl, Luca Fedeli
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This file is part of the WarpX automated test suite. It is used to test the
# injection of a laser pulse from an external binary file.
#
# - Generate an input binary file with a gaussian laser pulse.
# - Run the WarpX simulation for time T, when the pulse is fully injected
# - Compute the theory for laser envelope at time T
# - Compare theory and simulation, for both envelope and central frequency

import yt ; yt.funcs.mylog.setLevel(50)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import hilbert
import glob
import os

#Maximum acceptable error for this test
relative_error_threshold = 0.065

#Physical parameters
um = 1.e-6
fs = 1.e-15
c = 299792458

#Parameters of the gaussian beam
wavelength = 1.*um
w0 = 6.*um
tt = 10.*fs
x_c = 0.*um
t_c = 20.*fs
foc_dist = 10*um
E_max = 1e12
rot_angle = -np.pi/4.0

#Parameters of the tx grid
x_l = -12.0*um
x_r = 12.0*um
x_points = 480
t_l = 0.0*fs
t_r = 40.0*fs
t_points = 400
tcoords = np.linspace(t_l, t_r, t_points)
xcoords = np.linspace(x_l, x_r, x_points)

def gauss(T,X,Y,opt):
    """Compute the electric field for a Gaussian laser pulse.
       This is used to write the binary input file.
    """

    k0 = 2.0*np.pi/wavelength
    inv_tau2 = 1./tt/tt
    osc_phase = k0*c*(T-t_c)

    diff_factor = 1.0 + 1.0j* foc_dist * 2/(k0*w0*w0)
    inv_w_2 = 1.0/(w0*w0*diff_factor)

    pre_fact = np.exp(1.0j * osc_phase)

    if opt == '3d':
        pre_fact = pre_fact/diff_factor
    else:
        pre_fact = pre_fact/np.sqrt(diff_factor)

    exp_arg = - (X*X + Y*Y)*inv_w_2 - inv_tau2 * (T-t_c)*(T-t_c)

    return np.real(pre_fact * np.exp(exp_arg))

# Function for the envelope
def gauss_env(T,XX,ZZ):
    '''Function to compute the theory for the envelope
    '''

    X = np.cos(rot_angle)*XX + np.sin(rot_angle)*ZZ
    Z = -np.sin(rot_angle)*XX + np.cos(rot_angle)*ZZ

    inv_tau2 = 1./tt/tt
    inv_w_2 = 1.0/(w0*w0)
    exp_arg = - (X*X)*inv_w_2 - inv_tau2 / c/c * (Z-T*c)*(Z-T*c)
    return E_max * np.real(np.exp(exp_arg))

def write_file(fname, x, y, t, E):
    """ For a given filename fname, space coordinates x and y, time coordinate t
    and field E, write a WarpX-compatible input binary file containing the
    profile of the laser pulse
    """

    with open(fname, 'wb') as file:
        flag_unif = 0
        file.write(flag_unif.to_bytes(1, byteorder='little'))
        file.write((len(t)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(x)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(y)).to_bytes(4, byteorder='little', signed=False))
        file.write(t.tobytes())
        file.write(x.tobytes())
        file.write(y.tobytes())
        file.write(E.tobytes())


def write_file_unf(fname, x, y, t, E):
    """ For a given filename fname, space coordinates x and y, time coordinate t
    and field E, write a WarpX-compatible input binary file containing the
    profile of the laser pulse. This function should be used in the case
    of a uniform spatio-temporal mesh
    """

    with open(fname, 'wb') as file:
        flag_unif = 1
        file.write(flag_unif.to_bytes(1, byteorder='little'))
        file.write((len(t)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(x)).to_bytes(4, byteorder='little', signed=False))
        file.write((len(y)).to_bytes(4, byteorder='little', signed=False))
        file.write(t[0].tobytes())
        file.write(t[-1].tobytes())
        file.write(x[0].tobytes())
        file.write(x[-1].tobytes())
        if len(y) == 1 :
            file.write(y[0].tobytes())
        else :
            file.write(y[0].tobytes())
            file.write(y[-1].tobytes())
        file.write(E.tobytes())

def create_gaussian_2d():
   T, X, Y = np.meshgrid(tcoords, xcoords, np.array([0.0]), indexing='ij')
   E_t = gauss(T,X,Y,'2d')
   write_file("gauss_2d.txye", xcoords, np.array([0.0]), tcoords, E_t)
   write_file_unf("gauss_2d_unf.txye", xcoords, np.array([0.0]), tcoords, E_t)


def do_analysis(fname, compname, steps):
    ds = yt.load(fname)

    dt = ds.current_time.to_value()/steps

    # Define 2D meshes
    x = np.linspace(
        ds.domain_left_edge[0],
        ds.domain_right_edge[0],
        ds.domain_dimensions[0]).v
    z = np.linspace(
        ds.domain_left_edge[ds.dimensionality-1],
        ds.domain_right_edge[ds.dimensionality-1],
        ds.domain_dimensions[ds.dimensionality-1]).v
    X, Z = np.meshgrid(x, z, sparse=False, indexing='ij')

    # Compute the theory for envelope
    env_theory = gauss_env(+t_c-ds.current_time.to_value(), X,Z)+gauss_env(-t_c+ds.current_time.to_value(), X,Z)

    # Read laser field in PIC simulation, and compute envelope
    all_data_level_0 = ds.covering_grid(level=0,left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    F_laser = all_data_level_0['boxlib', 'Ey'].v.squeeze()
    env = abs(hilbert(F_laser))
    extent = [ds.domain_left_edge[ds.dimensionality-1], ds.domain_right_edge[ds.dimensionality-1],
          ds.domain_left_edge[0], ds.domain_right_edge[0] ]

    # Plot results
    plt.figure(figsize=(8,6))
    plt.subplot(221)
    plt.title('PIC field')
    plt.imshow(F_laser, extent=extent)
    plt.colorbar()
    plt.subplot(222)
    plt.title('PIC envelope')
    plt.imshow(env, extent=extent)
    plt.colorbar()
    plt.subplot(223)
    plt.title('Theory envelope')
    plt.imshow(env_theory, extent=extent)
    plt.colorbar()
    plt.subplot(224)
    plt.title('Difference')
    plt.imshow(env-env_theory, extent=extent)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(compname, bbox_inches='tight')

    relative_error_env = np.sum(np.abs(env-env_theory)) / np.sum(np.abs(env))
    print("Relative error envelope: ", relative_error_env)
    assert(relative_error_env < relative_error_threshold)

    fft_F_laser = np.fft.fft2(F_laser)

    freq_rows = np.fft.fftfreq(F_laser.shape[0],dt)
    freq_cols = np.fft.fftfreq(F_laser.shape[1],dt)

    pos_max = np.unravel_index(np.abs(fft_F_laser).argmax(), fft_F_laser.shape)

    freq = np.sqrt((freq_rows[pos_max[0]])**2 +  (freq_cols[pos_max[1]]**2))
    exp_freq = c/wavelength

    relative_error_freq = np.abs(freq-exp_freq)/exp_freq
    print("Relative error frequency: ", relative_error_freq)
    assert(relative_error_freq < relative_error_threshold)



def launch_analysis(executable):
    create_gaussian_2d()
    os.system("./" + executable + " inputs.2d_test_txye")
    do_analysis("diags/plotfiles/plt00250/", "comp_unf.pdf", 250)
    os.system("sed 's/gauss_2d_unf.txye/gauss_2d.txye/g' inputs.2d_test_txye > inputs.2d_test_txye_non_unf")
    os.system("./" + executable + " inputs.2d_test_txye_non_unf")
    do_analysis("diags/plotfiles/plt00250/", "comp_non_unf.pdf", 250)


def main() :
    executables = glob.glob("main2d*")
    if len(executables) == 1 :
        launch_analysis(executables[0])
    else :
        assert(False)
    print('Passed')

if __name__ == "__main__":
    main()
