#include <PlasmaInjector.H>
#include <WarpXConst.H>
#include <cmath>

#include <iostream>

using namespace amrex;
using namespace PhysConst;
using namespace MathConst;

/* This momentum distribution corresponds to a 3D plasma wave.
   The fields are given by the analytical formulas:
   $$ \phi = \epsilon \,\frac{m_e c^2}{q_e}\sin(k_x x)\sin(k_y y)\sin(k_z z) \sin(\omega_p t) $$
   $$ E_x = -\epsilon \,\frac{m_e c^2 k_x}{q_e}\cos(k_x x)\sin(k_y y)\sin(k_z z)\sin( \omega_p t)$$
   $$ E_y = -\epsilon \,\frac{m_e c^2 k_y}{q_e}\sin(k_x x)\cos(k_y y)\sin(k_z z)\sin( \omega_p t)$$
   $$ E_z = -\epsilon \,\frac{m_e c^2 k_z}{q_e}\sin(k_x x)\sin(k_y y)\cos(k_z z)\sin( \omega_p t)$$
   $$ v_x/c = \epsilon \,\frac{m_e q c k_x}{m q_e \omega_p}\cos(k_x x)\sin(k_y y)\sin(k_z z)\cos(\omega_p t) $$
   $$ v_y/c = \epsilon \,\frac{m_e q c k_y}{m q_e \omega_p}\sin(k_x x)\cos(k_y y)\sin(k_z z)\cos(\omega_p t) $$
   $$ v_z/c = \epsilon \,\frac{m_e q c k_z}{m q_e \omega_p}\sin(k_x x)\sin(k_y y)\cos(k_z z)\cos(\omega_p t) $$
*/

void CustomMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
  // The first parameter represents the amplitude of the plasma wave
  // The second parameter is the plasma density
  // The other parameters represents the number of periods in x, y, z
  Real epsilon = params[0];
  Real n = params[1];
  // NB: This assumes a box from -20 microns to 20 microns
  Real kx = 2.*pi*params[2]/40.e-6;
  Real ky = 2.*pi*params[3]/40.e-6;
  Real kz = 2.*pi*params[4]/40.e-6;
  // Plasma frequency
  Real wp = std::sqrt((n*q_e*q_e)/(m_e*ep0));

  // Because we are only considering electron and positrons here,
  // the ratio m_e * q / m q_e is either 1 or -1, and is encoded
  // in the sign of epsilon, as given in the input file
  u[0] = epsilon * c*kx/wp * std::cos(kx*x) * std::sin(ky*y) * std::sin(kz*z);
  u[1] = epsilon * c*kx/wp * std::sin(kx*x) * std::cos(ky*y) * std::sin(kz*z);
  u[2] = epsilon * c*kx/wp * std::sin(kx*x) * std::sin(ky*y) * std::cos(kz*z);
}
