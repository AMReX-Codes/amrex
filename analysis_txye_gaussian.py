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
w = 6.*um
tt = 10.*fs
x_c = 0.*um
y_c = 0.*um
t_c = 20.*fs
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

def gauss2d(T,X):
    return np.exp(-np.power(X-x_c,2)/(w*w))*np.sin(2.0*np.pi*c*(T-t_c)/wavelength)* np.exp(-np.power(T-t_c,2)/(tt*tt))


def gauss3d(T,X,Y):
    return np.exp(-np.power(Y-y_c,2)/(w*w))*gauss2d(T,X)


def test_func(x, a, b):
    return a * np.sin(b * x)


def analyze_gaussian_2d():
    filename = sys.argv[1]
    data_set_end = yt.load(filename)

    sim_time = data_set_end.current_time.to_value()

    ray0 = data_set_end.ray((0*um,0*um,0), (15*um, 15*um,0))
    ray1 = data_set_end.ray((-um/np.sqrt(2),um/np.sqrt(2),0), (-um/np.sqrt(2)+15*um, um/np.sqrt(2)+15*um,0))

    xx0 = np.array(ray0["t"])*np.sqrt(2)*15*um
    EE0 = np.array(ray0["Ey"])/E_max
    xx1 = np.array(ray1["t"])*np.sqrt(2)*12*um
    EE1 = np.array(ray1["Ey"])/E_max

    expected0 = [-gauss2d(x/c , 0) for x in xx0]
    #expected1 = [gauss2d(x/c, 1*um*np.sqrt(2)) for x in xx1]

    plt.plot(xx0,expected0)
    plt.plot(xx0,EE0,'bo')

    plt.show()


def main() :
    analyze_gaussian_2d()


if __name__ == "__main__":
    main()

