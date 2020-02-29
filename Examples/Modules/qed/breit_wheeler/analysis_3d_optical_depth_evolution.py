#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# -*- coding: utf-8 -*-

import yt
import numpy as np
import sys
import math as m
import scipy.special as spe
import scipy.integrate as integ

# This script checks if the optical depth of photons undergoing the
# Breit Wheeler process behaves as expected. Four populations of photons
# are initialized with different momenta in a background EM field. The decrease
# of the optical depth (averaged for each population) is used to estimate the
# Breit Wheeler cross section, which is then compared with the theoretical
# formula. Relative discrepancy should be smaller than 2%
#
# References:
# 1) R. Duclous et al 2011 Plasma Phys. Control. Fusion 53 015009
# 2) A. Gonoskov et al. 2015 Phys. Rev. E 92, 023305
# 3) M. Lobet. PhD thesis "Effets radiatifs et d'electrodynamique
#    quantique dans l'interaction laser-matiere ultra-relativiste"
#    URL: https://tel.archives-ouvertes.fr/tel-01314224


# Tolerance
tol = 2e-2

# EM fields
E_f = np.array([-2433321316961438, 973328526784575, 1459992790176863])
B_f = np.array([2857142.85714286, 4285714.28571428, 8571428.57142857])

# Physical constants
electron_mass = 9.10938356e-31
elementary_charge = 1.6021766208e-19
speed_of_light = 299792458
reduced_plank = 1.054571800e-34
vacuum_permittivity =  8.854187817e-12
fine_structure_constant =  0.0072973525664
classical_elec_radius = (1./4./np.pi/vacuum_permittivity)*( elementary_charge**2 / (electron_mass * speed_of_light**2))
lambda_ref = 1.0e-6
field_reference = 2.0 * np.pi * electron_mass*speed_of_light * speed_of_light / (elementary_charge*lambda_ref)
schwinger_field_SI = electron_mass**2 * speed_of_light**3/(reduced_plank*elementary_charge)
schwinger_field_norm = electron_mass*speed_of_light*lambda_ref/(2.0*reduced_plank*m.pi)
#______________

# Initial parameters
spec_names = ["p1", "p2", "p3", "p4"]
 #Assumes average = 1 at t = 0 (ok for a large number of particles)
tau_begin_avg = np.array([1.0, 1.0, 1.0, 1.0])
mec = electron_mass*speed_of_light
p_begin = [
    np.array([2000.0,0,0])*mec,
    np.array([0.0,5000.0,0.0])*mec,
    np.array([0.0,0.0,10000.0])*mec,
    np.array([57735.02691896, 57735.02691896, 57735.02691896])*mec
]
#______________

def calc_chi_gamma(p, E, B):
    p = p / electron_mass / speed_of_light
    E = E / field_reference
    B = B * speed_of_light / field_reference
    gamma_phot = np.linalg.norm(p)
    c = p/gamma_phot
    loc_field = gamma_phot * np.linalg.norm( E - np.dot(c,E)*c + np.cross(c,B))
    return loc_field/schwinger_field_norm

#Auxiliary functions
def X(chi_phot, chi_ele):
    if (chi_phot > chi_ele and chi_ele != 0):
        return np.power(chi_phot/(chi_ele*(chi_phot-chi_ele)), 2./3.)
    else:
        return 1.0e30

def T(chi_phot):
    coeff = 1./(np.pi * np.sqrt(3.) * chi_phot * chi_phot)
    inner = lambda x : integ.quad(lambda s: np.sqrt(s)*spe.kv(1./3., 2./3. * s**(3./2.)), x, np.inf)[0]
    return integ.quad(lambda chi_ele:
                      coeff*(inner(X(chi_phot, chi_ele)) -
                      (2.0 - chi_phot*np.power(X(chi_phot, chi_ele), 3./2.))*spe.kv(2./3., 2./3. *X(chi_phot, chi_ele)**(3./2.)) )
                      , 0, chi_phot)[0]
#__________________

# Breit-Wheeler total cross section
def dNBW_dt(chi_phot, energy_phot):
    energy_phot = energy_phot/electron_mass/speed_of_light/speed_of_light
    return ((electron_mass*(speed_of_light)**2)*fine_structure_constant/reduced_plank)*(chi_phot/energy_phot)*T(chi_phot)
#__________________

def check():
    filename_end = sys.argv[1]
    data_set_end = yt.load(filename_end)

    sim_time = data_set_end.current_time.to_value()

    all_data_end = data_set_end.all_data()

    tau_end_avg = np.array([
       np.average(all_data_end[name, 'particle_optical_depth_BW'])
       for name in spec_names])

    dNBW_dt_sim = (tau_begin_avg - tau_end_avg)/sim_time

    dNBW_dt_theo = [
        dNBW_dt(calc_chi_gamma(p, E_f, B_f), np.linalg.norm(p*speed_of_light))
        for p in p_begin
    ]

    discrepancy = np.abs(dNBW_dt_sim-dNBW_dt_theo)/dNBW_dt_theo

    assert(np.all(discrepancy < tol))

def main():
    check()

if __name__ == "__main__":
    main()
