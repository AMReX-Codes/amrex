# coding: utf-8
"""Generate and test numerical routines for a 2d Chemotaxis PDE.

The "minimal" Keller-Segel chemotaxis model is

  u_t = ∇·( D ∇u + Χ ∇v )
  v_t = ∇·( ∇v ) + u - v

where u is the cell (organism) density and v is the concentration of
the chemical signal.

For the u_t equation, D ∇u is the cell/organism "motility", and Χ ∇v
is the "chemotactic sensitivity".

"""

import numpy as np
import sympy

import codegen.compile
import codegen.symbolic
from codegen.polyquad import Kernel
from codegen.symbolic import Grad as grad

from numpy import pi
from sympy.abc import x, y

from ctypes import c_int, c_double
from numpy.ctypeslib import ndpointer

c_int_2 = 2*c_int

#### code generation

# define unknowns and constants

u = sympy.Function('u')
v = sympy.Function('v')

n     = 1
diff  = sympy.var('diff')
chi   = sympy.var('chi')
alpha = sympy.var('alpha')
gamma = sympy.var('gamma')
phi   = sympy.var('phi')

# cell motility (flux)

operator = diff * grad(u(x, y))    # minimal model
operator = operator * u(x, y)**n   # non-linear diffusion

motility = Kernel()
motility.gen_flux(operator)

# chemotactic sensitivity (flux)

operator = chi * u(x,y) * grad(v(x,y))            # minimal model
operator = operator * (1 - u(x, y)/gamma)         # volume filling
operator = operator / (1 + alpha * v(x, y))**2    # receptor binding

sensitivity = Kernel()
sensitivity.gen_flux(operator)

# signal diffusion (flux)

operator = grad(v(x,y))                     # minimal model

diffusion = Kernel()
diffusion.gen_flux(operator)

# signal production (source)

operator = u(x,y)                           # minimal model
operator = operator / (1 + phi * u(x, y))   # non-linear signal kinetics

production = Kernel()
production.gen_source(operator)

# write to kernels.f90

with open("kernels.tmpl.f90", 'r') as f:
    tmpl = f.read()

with open("kernels.f90", 'w') as f:
    f.write(tmpl.format(
            motility_kernel=motility.flux(),
            sensitivity_kernel=sensitivity.flux(),
            diffusion_kernel=diffusion.flux(),
            production_kernel=production.source(),
            ))

# compile and set prototypes

module = codegen.compile.compile("kernels")

module.cell_motility.argtypes = [ ndpointer(ndim=2, dtype='f8'), 
                                  ndpointer(ndim=2, dtype='f8'), 
                                  c_int_2, c_int_2, c_int, c_double, c_double ]

module.chemotactic_sensitivity.argtypes = [ ndpointer(ndim=2, dtype='f8'), 
                                            ndpointer(ndim=2, dtype='f8'), 
                                            ndpointer(ndim=2, dtype='f8'),
                                            c_int_2, c_int_2, c_int, c_double, c_double, c_double, c_double ]

module.signal_production.argtypes = [ ndpointer(ndim=2, dtype='f8'), 
                                      ndpointer(ndim=2, dtype='f8'),
                                      ndpointer(ndim=2, dtype='f8'),
                                      c_int_2, c_int_2, c_int, c_double, c_double ]

raise SystemExit

#### validation of generated code

# define test functions

def u(x, y):
    return 1.0 + 0.1 * sympy.cos(pi*x) * sympy.cos(2*pi*y)

def v(x, y):
    return sympy.cos(2*pi*x)

def f_motility(x, y):
    return symbolic.grad(u(x, y))

def f_sensitivity(x, y):
    return u(x,y)*symbolic.grad(v(x, y))

def f_production(x, y):
    return u(x, y) / (1 + u(x, y))

# compute cell averaged u and v, and exact fluxes for the motility and
# sensitivity terms

ng   = 2
nx   = 16
dx   = 2.0/nx
grid = np.mgrid[-1.0:1.0:(nx+1)*1j, -1.0:1.0:(nx+1)*1j]

ubar = symbolic.mesh_volume_integrals_2d(u, grid) / dx**2
vbar = symbolic.mesh_volume_integrals_2d(v, grid) / dx**2

# motility_exact    = symbolic.mesh_surface_integrals_2d(f_motility, grid)
# sensitivity_exact = symbolic.mesh_surface_integrals_2d(f_sensitivity, grid)
production_exact = symbolic.mesh_volume_integrals_2d(f_production, grid)

from pyboxlib import *
pybl.open()

la = layout()
la.create(boxes=[ ( (1, 1), (nx, nx) ) ])
 
U = multifab() 
U.create(la, components=2, ghost_cells=ng)

fab = U.fab(1)
fab[1:nx+1, 1:nx+1, 0] = ubar
fab[1:nx+1, 1:nx+1, 1] = vbar

U.fill_boundary()

F = np.zeros((nx, nx), order='F')

lo, hi = c_int_2(), c_int_2()
lo[0], lo[1] = 1, 1
hi[0], hi[1] = nx, nx

#module.cell_motility(F, fab.array[:, :, 0], lo, hi, 2, 2.0/nx, 1.0)
#module.chemotactic_sensitivity(F, fab.array[:, :, 0], fab.array[:, :, 1], lo, hi, ng, 2.0/nx, 1.0)
module.signal_production(F, fab.array[:, :, 0], lo, hi, ng, 2.0/nx, 1.0)

import matplotlib.pylab as plt

plt.figure(1)
plt.imshow(production_exact.transpose())
plt.colorbar()

plt.figure(3)
plt.imshow(F.transpose())
plt.colorbar()

#print np.log10(abs(motility_exact - F).max())
#print np.log10(abs(sensitivity_exact - F).max())
print np.log10(abs(production_exact - F).max())

plt.show()


