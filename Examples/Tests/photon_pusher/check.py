#! /usr/bin/env python
import yt
import numpy as np
import sys

filename = sys.argv[1]
data_set_end = yt.load(filename)
print(filename)

sim_time = data_set_end.current_time.to_value()

c = 299792458.0
m_e = 9.10938291e-31

#expected positions
ll = sim_time*c
xp = np.array([sim_time*c, 0.0, 0.0])
xn = np.array([-sim_time*c, 0.0, 0.0])
yp = np.array([0.0, sim_time*c, 0.0])
yn = np.array([0.0, -sim_time*c, 0.0])
zp = np.array([0.0, 0.0, sim_time*c])
zn = np.array([0.0, 0.0, -sim_time*c])
dp = np.array([sim_time*c/np.sqrt(3.0), sim_time*c/np.sqrt(3.0), sim_time*c/np.sqrt(3.0)])
dn = np.array([-sim_time*c/np.sqrt(3.0), -sim_time*c/np.sqrt(3.0), -sim_time*c/np.sqrt(3.0)])

#expected positions list
answ_pos = [xp, xn, yp, yn, zp, zn, dp, dn, xp, xn, yp, yn, zp, zn, dp, dn]

#espected momenta
pp1 =  m_e * c * 1
pp10 = m_e * c * 10
mxp1 = np.array([pp1, 0.0, 0.0])
mxn1 = np.array([-pp1, 0.0, 0.0])
myp1 = np.array([0.0, pp1, 0.0])
myn1 = np.array([0.0, -pp1, 0.0])
mzp1 = np.array([0.0, 0.0, pp1])
mzn1 = np.array([0.0, 0.0, -pp1])
mdp1 = np.array([pp1, pp1, pp1])
mdn1 = np.array([-pp1, -pp1, -pp1])
mxp10 = np.array([pp10, 0.0, 0.0])
mxn10 = np.array([-pp10, 0.0, 0.0])
myp10 = np.array([0.0, pp10, 0.0])
myn10 = np.array([0.0, -pp10, 0.0])
mzp10 = np.array([0.0, 0.0, pp10])
mzn10 = np.array([0.0, 0.0, -pp10])
mdp10 = np.array([pp10, pp10, pp10])
mdn10 = np.array([-pp10,-pp10, -pp10])

#expected momenta list
answ_mom = [mxp1, mxn1, myp1, myn1, mzp1, mzn1, mdp1, mdn1,
mxp10, mxn10, myp10, myn10, mzp10, mzn10, mdp10, mdn10]


#species names
spec_names = ["p_xp_1", "p_xn_1", "p_yp_1", "p_yn_1",
"p_zp_1", "p_zn_1","p_dp_1", "p_dn_1",
"p_xp_10", "p_xn_10", "p_yp_10", "p_yn_10",
"p_zp_10", "p_zn_10", "p_dp_10", "p_dn_10"]

#simulation results
all_data =  data_set_end.all_data()
res_pos = [np.array([
all_data[sp, 'particle_position_x'].v[0],
all_data[sp, 'particle_position_y'].v[0],
all_data[sp, 'particle_position_z'].v[0]])
for sp in spec_names]
res_mom = [np.array([
all_data[sp, 'particle_momentum_x'].v[0],
all_data[sp, 'particle_momentum_y'].v[0],
all_data[sp, 'particle_momentum_z'].v[0]])
for sp in spec_names]

#check discrepancies
tol_pos = 1.0e-14;
tol_mom = 0.0; #momentum should be conserved exactly
disc_pos = [np.linalg.norm(a-b)/np.linalg.norm(b)
for a,b in zip(res_pos, answ_pos)]
disc_mom = [np.linalg.norm(a-b)/np.linalg.norm(b)
for a,b in zip(res_mom, answ_mom)]

assert ((max(disc_pos) <= tol_pos) and (max(disc_mom) <= tol_mom))
