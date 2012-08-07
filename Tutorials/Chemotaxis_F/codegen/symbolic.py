# coding: utf-8
"""Symbolic helper functions for testing and constructing FV schemes."""

import numpy as np
import sympy

from sympy.abc import x, y, z

dim = 2

def grad(f, dim=dim, vars=[x,y,z]):
    """Compute the gradient of f."""
    vars = vars[:dim]
    grad = sympy.zeros(len(vars), 1)
    for i in range(len(vars)):
        grad[i, 0] = f.diff(vars[i])
    return grad


class Grad(sympy.Basic):

    def __new__(cls, f, vars=[x, y]):
        f = sympy.sympify(f)
        r = cls.canonicalize(f, vars)
        if r is not None:
            return r
        return Basic.__new__(cls, f, vars)

    @classmethod
    def canonicalize(cls, f, vars):
        grad = sympy.zeros(len(vars), 1)
        for i in range(len(vars)):
            grad[i, 0] = sympy.diff(f, vars[i])
        return grad

    def tostr(self):
        return "grad(%s, %s)" % (self.args[0], self.args[1])



# def div(f, dim=dim, vars=[x,y,z]):
#     """Compute the divergence of f."""
#     vars = vars[:dim]
#     div = 0
#     for i in range(len(vars)):
#         div += f[i].diff(vars[i])
#     return div


_si_cache = {}

def surface_integral_2d(f, xs=[0.0, 1.0], ys=[0.0, 1.0], numeric=True):
    """Compute the surface integral of f·dA over a box.""" 

    from sympy import integrate
    from sympy.mpmath import quad
    from sympy.abc import x, y

    f = f(x, y)

    if numeric:
        i1 = quad( f[0].subs(x, xs[1]), (y, ys[0], ys[1]))
        i2 = quad( f[1].subs(y, ys[1]), (x, xs[0], xs[1]))
        i3 = quad(-f[0].subs(x, xs[0]), (y, ys[0], ys[1]))
        i4 = quad(-f[1].subs(y, ys[0]), (x, xs[0], xs[1]))

    else:

        if f[0] not in _si_cache:
            _si_cache[f[0]] = integrate(f[0], y)

        if f[1] not in _si_cache:
            _si_cache[f[1]] = integrate(f[1], x)

        F0 = _si_cache[f[0]]
        F1 = _si_cache[f[1]]

        # i1 = integrate( f[0].subs(x, xs[1]), (y, ys[0], ys[1]))
        i1 = F0.subs(x, xs[1])
        i1 = i1.subs(y, ys[1]) - i1.subs(y, ys[0])

        # i2 = integrate( f[1].subs(y, ys[1]), (x, xs[0], xs[1]))
        i2 = F1.subs(y, ys[1])
        i2 = i2.subs(x, xs[1]) - i2.subs(x, xs[0])

        # i3 = integrate(-f[0].subs(x, xs[0]), (y, ys[0], ys[1]))
        i3 = -F0.subs(x, xs[0])
        i3 = i3.subs(y, ys[1]) - i3.subs(y, ys[0])

        # i4 = integrate(-f[1].subs(y, ys[0]), (x, xs[0], xs[1]))
        i4 = -F1.subs(y, ys[0])
        i4 = i4.subs(x, xs[1]) - i4.subs(x, xs[0])

    return i1 + i2 + i3 + i4


def volume_integral_2d(f, xs=[0.0, 1.0], ys=[0.0, 1.0]):
    """Compute the volume integral of f dV over a box."""

    from sympy import integrate
    from sympy.mpmath import quad
    from sympy.abc import x, y

    # f = f(x, y)
    # i1 = integrate(f,  (y, ys[0], ys[1]))
    # i2 = integrate(i1, (x, xs[0], xs[1]))

    g = lambda x, y: f(x, y).evalf()

    i2 = quad(g, [ xs[0], xs[1] ], [ ys[0], ys[1] ] )
      
    return i2


def mesh_surface_integrals_2d(f, mesh, verbose=False):
    """Compute the surface integral of f over each grid box in the mesh."""

    dim, nx, ny = mesh.shape
    nx, ny = nx-1, ny-1

    flux = np.zeros((nx, ny))

    for i in range(nx):
        for j in range(ny):
            if verbose:
                print 'surface integral', i, j
            flux[i, j] = surface_integral_2d(f, 
                              xs=[ mesh[0, i, j], mesh[0, i+1, j] ],
                              ys=[ mesh[1, i, j], mesh[1, i, j+1] ])

    return flux

def mesh_volume_integrals_2d(f, mesh):
    """Compute the volume integral of f over each grid box in the mesh."""

    dim, nx, ny = mesh.shape
    nx, ny = nx-1, ny-1

    avg = np.zeros((nx, ny))

    for i in range(nx):
        for j in range(ny):
            avg[i, j] = volume_integral_2d(f, 
                              xs=[ mesh[0, i, j], mesh[0, i+1, j] ],
                              ys=[ mesh[1, i, j], mesh[1, i, j+1] ])

    return avg
  

