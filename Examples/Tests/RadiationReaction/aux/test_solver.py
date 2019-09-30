# This file contains few simple tests of solver.py.
# It initializes an electron with normalized momentum in different
# directions, propagating in a static magnetic field having normalized
# amplitude equal to 5975 (lambda = 1 um)
# M. Vranic et al. Computer Physics Communications 204 (2016):
# 141-151. https://arxiv.org/pdf/1502.02432.pdf ) states that
# if the initial momentum is perpendicular to the magnetic field,
# after 16 giration times the final total energy should be ~ 6 mc2.
# If the initial momentum is parallel to the field, it should not change.
# This module tests these predictions respectively with a tolerance of
# 5% and 1e-4%

import solver
import numpy as np

# TESTS

def case1():
    rel_tol =  0.05
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [pval,0.,0.]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([0,0,Bval])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = 6.0
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case2():
    rel_tol =  0.05
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [0.0,pval,0.]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([0,0,Bval])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = 6.0
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case3():
    rel_tol =  1e-6
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [0.0,0.,pval]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([0,0,Bval])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = np.sqrt(1.0 + pval*pval)
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case4():
    rel_tol =  0.05
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [pval,0.,0.]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([0,Bval,0])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = 6.0
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case5():
    rel_tol =  0.05
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [0.0,0.,pval]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([0,Bval,0])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = 6.0
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case6():
    rel_tol =  1e-6
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [0.0,pval,0.0]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([0,Bval,0])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = np.sqrt(1.0 + pval*pval)
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case7():
    rel_tol =  0.05
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [0.0,pval,0.]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([Bval,0.0,0])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = 6.0
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case8():
    rel_tol =  0.05
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [0.0,0.,pval]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([Bval,0.0,0])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = 6.0
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

def case9():
    rel_tol =  1e-6
    init_pos = [0.,0.,0.]
    pval = 10.0
    Bval = 5975
    init_mom = [pval,0.0, 0.0]
    init_t = 0.0
    Tgir =  2.0*np.pi*np.sqrt(1. + pval*pval)/Bval
    end_t = 16.0*Tgir
    steps = 20000
    E = lambda pos, t : np.array([0,0,0])
    B = lambda pos, t : np.array([Bval,0.0,0])
    charge = -1
    ll = 1.0e-6
    sol = solver.solve(init_pos, init_mom, init_t, end_t, steps, E, B, charge, ll)

    exp_res = np.sqrt(1.0 + pval*pval)
    res =  np.sqrt(1.0 + sol[-1,3]* sol[-1,3] + sol[-1,4]* sol[-1,4] + sol[-1,5]* sol[-1,5])

    assert(np.abs(res-exp_res)/res < rel_tol)

# ____________________________________________________________________________

def do_test():
    case1()
    case2()
    case3()
    case4()
    case5()
    case6()
    case7()
    case8()
    case9()
