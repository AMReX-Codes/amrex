#!/usr/bin/python3

import sys
import numpy as np

#Physical parameters
um = 1.e-6
fs = 1.e-15
c = 299792458

#Parameters of the gaussian beam
wavelength = 1.*um
w = 6.*um
tt = 10.*fs
x_c = 0.*um
y_c = 0.*um
t_c = 20.*fs
foc_dist = 10.0*um

#Parameters of the txy grid
x_l = -12.0*um
x_r = 12.0*um
x_points = 480
y_l = -12.0*um
y_r = 12.0*um
y_points = 480
t_l = 0.0*fs
t_r = 40.0*fs
t_points = 200
tcoords = np.linspace(t_l, t_r, t_points)
xcoords = np.linspace(x_l, x_r, x_points)
ycoords = np.linspace(y_l, y_r, y_points)

def gauss(T,X,Y,opt):
    k0 = 2.0*np.pi/wavelength
    inv_tau2 = 1./tt/tt
    osc_phase = k0*c*(T-t_c)

    diff_factor = 1.0 + 1.0j* foc_dist * 2/(k0*w*w)
    inv_w_2 = 1.0/(w*w*diff_factor)

    pre_fact = 1.0*np.exp(1.0j * osc_phase)

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

def create_file_gnu(x, t, E):
    with open("data.dat", 'w') as file:
        print("gnu gnu")
        print(len(t), len(x))
        print(E.shape)
        print("gnu gnu")
        for i in range(len(t)):
            for j in range(len(x)):
                k = 0
                file.write("{} {} {} \n".format(t[i], x[j], E[i,j,k]))


def create_gaussian_2d():
   T, X, Y = np.meshgrid(tcoords, xcoords, np.array([0.0]), indexing='ij')
   E_t = gauss(T,X,Y,'2d')
   write_file("gauss_2d.txye", xcoords, np.array([0.0]), tcoords, E_t)
   write_file_unf("gauss_2d_unf.txye", xcoords, np.array([0.0]), tcoords, E_t)
   create_file_gnu(xcoords, tcoords, E_t)


def create_gaussian_3d():
   T, X, Y = np.meshgrid(tcoords, xcoords, ycoords, indexing='ij')
   E_t = gauss(T,X,Y,'3d')
   write_file("gauss_3d.txye", xcoords, ycoords, tcoords, E_t)
   write_file_unf("gauss_3d_unf.txye", xcoords, ycoords, tcoords, E_t)


def main() :
    create_gaussian_2d()


if __name__ == "__main__":
    main()

