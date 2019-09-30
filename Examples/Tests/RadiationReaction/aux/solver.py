# An odepack-based solver to integrate Reduced Landau-Lifshitz equations
# See
# M. Vranic et al. Computer Physics Communications 204 (2016):
# 141-151. https://arxiv.org/pdf/1502.02432.pdf )
# and
# M Tamburini et al 2010 New J. Phys. 12 (2010):
# 123005 https://iopscience.iop.org/article/10.1088/1367-2630/12/12/123005/meta

import numpy as np
import scipy.integrate as scinteg

classical_electron_radius = 2.81794e-15

def RLL (yy, t, E, B, charge, ll):
    pos = yy[0:3]
    mom = yy[3:6]

    gamma2 = 1. + np.dot(mom, mom)
    inv = 1./np.sqrt(gamma2)

    v = mom*inv

    EE = E(pos, t)
    BB = B(pos, t)

    f_l = charge*(EE + np.cross(v, BB))

    coeff = - (4./3.)*np.pi * classical_electron_radius/ll

    vdotE = np.dot(v, EE)
    t1 = np.cross(f_l, BB) - vdotE*EE
    t2 = gamma2*(np.dot(f_l, f_l) - vdotE*vdotE)*v
    f_rr = coeff*(t1 + t2)

    yydot_pos = v
    yydot_mom = f_l + f_rr

    return np.concatenate((yydot_pos, yydot_mom))


def solve(init_pos, init_mom, init_time, end_time, nsteps, E, B, charge, ll):
    Y0 = np.concatenate((init_pos, init_mom))

    t = np.linspace(init_time, end_time, num=nsteps)

    func = lambda yy, tt : RLL(yy, tt, E, B, charge, ll)
    return scinteg.odeint(func, Y0, t)
