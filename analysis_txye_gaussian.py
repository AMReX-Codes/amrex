#!/usr/bin/python3

import yt
import sys
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

#Physical parameters
um = 1.e-6
fs = 1.e-15
c = 299792458

#Parameters of the gaussian beam
wavelength = 1.*um
w0 = 6.*um
tt = 10.*fs
x_c = 0.*um
y_c = 0.*um
t_c = 20.*fs
foc_dist = 10*um
E_max = 1e12

#Parameters of the txy grid
x_l = -5.0*um
x_r = 5.0*um
x_points = 200
y_l = -5.0*um
y_r = 5.0*um
y_points = 200
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

    diff_factor = 1.0 + 1.0j* foc_dist * 2/(k0*w0*w0)
    inv_w_2 = 1.0/(w0*w0*diff_factor)

    pre_fact = np.exp(1.0j * osc_phase)

    if opt == '3d':
        pre_fact = pre_fact/diff_factor
    else:
        pre_fact = pre_fact/np.sqrt(diff_factor)

    exp_arg = - (X*X + Y*Y)*inv_w_2 - inv_tau2 * (T-t_c)*(T-t_c)

    return np.real(pre_fact * np.exp(exp_arg))

def my_gauss(T,X,Y,opt):
    omega = 2.0*np.pi*c/wavelength
    Zr = omega * w0**2/2.
    w  = np.sqrt(1./(1.+(foc_dist/Zr)**2))
    invWaist2 = (w/w0)**2
    inv_tau2 = 1./tt/tt
    coeff = -omega * foc_dist * w**2 / (2.*Zr**2)
    space = np.exp( -invWaist2*((X-x_c)**2 + (Y-y_c)**2 ))
    if opt is '2d' :
        space = np.sqrt(w)*space
    else :
        space = w*space
    phase = coeff * ((X-x_c)**2 + (Y-y_c)**2)
    osc_phase =  (2.0*np.pi*c/wavelength)*(T-t_c)
    return np.real(space*np.exp( 1.0j*(osc_phase + phase) - inv_tau2 * (T-t_c)*(T-t_c)))

def my_gauss_fit(T,X,Y,opt, foc, wav, wa):
    omega = 2.0*np.pi*c/wav
    Zr = omega * wa**2/2.
    w  = np.sqrt(1./(1.+(foc/Zr)**2))
    invWaist2 = (w/wa)**2
    inv_tau2 = 1./tt/tt
    coeff = -omega * foc * w**2 / (2.*Zr**2)
    space = np.exp( -invWaist2*((X-x_c)**2 + (Y-y_c)**2 ))
    if opt is '2d' :
        space = np.sqrt(w)*space
    else :
        space = w*space
    phase = coeff * ((X-x_c)**2 + (Y-y_c)**2)
    osc_phase =  (2.0*np.pi*c/wav)*(T-t_c)
    return np.real(space*np.exp( 1.0j*(osc_phase + phase) - inv_tau2 * (T-t_c)*(T-t_c)))

def test_func(x, a, b):
    return a * np.sin(b * x)


def analyze_gaussian_2d():
    filename = sys.argv[1]
    data_set_end = yt.load(filename)

    sim_time = data_set_end.current_time.to_value()

    print(sim_time)
    dt = sim_time/507

    ray0 = data_set_end.ray((0*um,0*um,0), (15*um, 15*um,0))
    ray1 = data_set_end.ray((-um/np.sqrt(2),um/np.sqrt(2),0), (-um/np.sqrt(2)+15*um, um/np.sqrt(2)+15*um,0))

    xx0 = np.array(ray0["t"])*np.sqrt(2)*15*um
    EE0 = np.array(ray0["Ey"])/E_max
    xx1 = np.array(ray1["t"])*np.sqrt(2)*15*um
    EE1 = np.array(ray1["Ey"])/E_max


    expected0 = [-my_gauss((sim_time-0.5*dt)-x/c , 0, 0, '2d') for x in xx0]
    expected1 = [-my_gauss((sim_time)-x/c, 1*um*np.sqrt(2),0,'2d') for x in xx1]

    #plt.plot(xx0,expected0, 'ro')
    plt.plot(xx0,EE0,'bo')
    plt.plot(xx0,expected0-EE0, 'ro')
    #plt.plot(xx1,EE1,'bo')

    '''
    func = lambda x, foc, wav : my_gauss_fit((sim_time + dt*1)-x/c,0,0,opt, foc, wav, 1.0e-6)
    popt, pcov = opt.curve_fit(func, xx0, EE0,  bounds=([5.e-6, 1.0e-6], [20.e-6, 10.0e-6]))
    print(popt)
    print(pcov)
    plt.plot(xx0, func(xx0, *popt), 'ro')
    '''
    plt.show()


def main() :
    analyze_gaussian_2d()


if __name__ == "__main__":
    main()

