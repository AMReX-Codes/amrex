#!/usr/bin/python3

#COMMENT TO ADD
#COMMENT TO ADD
#COMMENT TO ADD
#COMMENT TO ADD

import sys
import numpy as np
import matplotlib.pyplot as plt
import yt

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

def write_file(fname, x, y, t, E):
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
        file.write(y[0].tobytes())
        file.write(y[-1].tobytes())
        file.write(E.tobytes())


def create_gaussian_2d():
   T, X, Y = np.meshgrid(tcoords, xcoords, np.array([0.0]), indexing='ij')
   E_t = gauss(T,X,Y,'2d')
   write_file("gauss_2d.txye", xcoords, np.array([0.0]), tcoords, E_t)
   write_file_unf("gauss_2d_unf.txye", xcoords, np.array([0.0]), tcoords, E_t)

def do_analysis(fname):
    data_set_end = yt.load(fname)
    sim_time = data_set_end.current_time.to_value()
    ray0 = data_set_end.ray((0*um,0*um,0), (15*um, 15*um,0))

    xx0 = np.array(ray0["t"])*np.sqrt(2)*15*um
    EE0 = np.array(ray0["Ey"])/E_max

    expected0 = [-gauss((sim_time)-x/c , 0, 0, '2d') for x in xx0]

    #DEBUG
    plt.plot(xx0,EE0,'bo')
    plt.plot(xx0,expected0,'ro')
    plt.savefig('graph.png')
    #__DEBUG__

    return True

def main() :
    arg = sys.argv[1]
    if arg == "--generate_txye_files":
        create_gaussian_2d()
    elif arg == "--check":
        do_analysis(sys.argv[2])


if __name__ == "__main__":
    main()

