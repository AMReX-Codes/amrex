"""Routines to build kernels by approximating operators with Gaussian
quadrature + polynomial reconstruction."""

import sympy
from sympy.abc import x, y


h = sympy.var('dx')


###############################################################################

def polynomial_interpolator(x, y):
  """Build a symbolic polynomial that interpolates the points (x_i, y_i).

  The returned polynomial is a function of the SymPy variable x.
  """

  xi = sympy.var('x')

  # build polynomial interpolant in Lagrange form
  k = len(x)
  sum_i = 0
  for i in range(k):

    ns = range(0,k)
    ns.remove(i)

    prod_n = 1
    for n in ns:
      prod_n = prod_n * (xi - x[n])

    prod_m = 1
    for m in ns:
      prod_m = prod_m * (x[i] - x[m])

    sum_i = sum_i + prod_n / prod_m * y[i]

  return sum_i


###############################################################################

def primitive_polynomial_interpolator(x, y):
  """Build a symbolic polynomial that approximates the primitive
  function f such that f(x_i) = sum_j y_j * (x_{j+1} - x_{j}).

  The returned polynomial is a function of the SymPy variable 'x'.
  """

  Y = [0]
  for i in range(0, len(y)):
    Y.append(Y[-1] + (x[i+1] - x[i]) * y[i])

  return polynomial_interpolator(x, Y)


###############################################################################

def _pt(a, b, x):
  """Map x in [-1, 1] to x in [a, b]."""

  half = sympy.sympify('1/2')

  w = half * (b - a)
  c = half * (a + b)

  return w * x + c


def reconstruction_coefficients(xi, ndiff=1, k=5):
  """Compute the polynomial reconstruction coefficients for the
  reconstruction points in *xi*."""

  (x, dx) = sympy.var('x dx')

  # cell boundaries
  xs = []
  for j in range(-k/2+1, k/2+2):
    xs.append(j*dx)

  # cell averages
  fs = []
  for j in range(-k/2+1, k/2+1):
    fs.append(sympy.var('f[i%+d]' % j))

  # compute reconstruction coefficients
  c = {}
  p = primitive_polynomial_interpolator(xs, fs).diff(x, ndiff)
  p = p.expand()

  for l in range(len(xi)):
      z = _pt(xs[k/2], xs[k/2+1], xi[l])

      for j in range(k):
          c[l, -k/2+j+1] = p.subs(x, z).coeff(fs[j])
          if c[l, -k/2+j+1] is None:
              c[l, -k/2+j+1] = 0.0

  return c


###############################################################################

def make_unknown(name):
    symbol = lambda x: sympy.Symbol(x, real=True)

    phi = {}
    for i in range(-5, 6):
        for j in range(-5, 6):
            phi[i, j] = symbol('%s(i%+d,j%+d)' % (name, i, j))
    return phi


def set_vec(vec, values):
    vals = []
    for v in values:
        v = sympy.sympify(v)
        v = sympy.fcode((v.evalf())).strip()
        if v == '0':
            v = '0.0d0'
        vals.append(v)
    return vec + ' = [ ' + ', '.join(vals) + ' ]'


def set_coeffs(c, xi, stencil):
    coeffs = [ c[xi, s] for s in stencil ]
    return set_vec('cof', coeffs)


def name_and_unknown(expr):
    l = []
    if isinstance(expr, sympy.Derivative):
        f = str(expr.args[0].__class__)
        v = str(expr.args[1])
        l.append((f + '_' + v, f, 1))
    elif isinstance(expr, sympy.Function):
        v = str(expr.__class__)
        l.append((v, v, 0))
    else:
        for arg in expr.args:
            l.extend(name_and_unknown(arg))
    return l


def expr_and_name(expr):
    l = []
    if isinstance(expr, sympy.Derivative):
        f = str(expr.args[0].__class__)
        v = str(expr.args[1])
        l.append((expr, f + '_' + v))
    elif isinstance(expr, sympy.Function):
        v = str(expr.__class__)
        l.append((expr, v))
    else:
        for arg in expr.args:
            l.extend(expr_and_name(arg))
    return l


def to_string_op(op, suffix):

    # first pass, only derivatives
    for expr, name in expr_and_name(op):
        if isinstance(expr, sympy.Derivative):
            op = op.subs(expr, name + suffix)

    # second pass, everything else
    for expr, name in expr_and_name(op):
        op = op.subs(expr, name + suffix)

    return str(op)



###############################################################################
# code generation routines

class Kernel(object):

    def __init__(self, k=5):

        self.k = 5
        self.stencil  = range(-k/2+1, k/2+1) 
        self.quad_pts = [ -sympy.sqrt(3)/5, 0, sympy.sqrt(3)/5 ]
        self.quad_wts = [ sympy.sympify('5/9'), sympy.sympify('8/9'), sympy.sympify('5/9') ]

        self.fsrc = []
        self.ssrc = []


    def set_quads(self):

        self.fsrc.append(set_vec('qwts', self.quad_wts))
        self.ssrc.append(set_vec('qwts', self.quad_wts))


    def reconstruct_edge(self, dest, unknown, edge, normal_derivative=0):
        """Generate code to reconstruct *unknown* along edge *edge*."""

        if (dest, edge) in self.rcache:
            return
        else:
            self.rcache[dest, edge] = True

        src = []
        phi = make_unknown(unknown)
        stencil = self.stencil

        if edge in [ 'right', 'top' ]:
            c = reconstruction_coefficients([ +1 ], k=self.k, ndiff=normal_derivative+1)
        else: 
            c = reconstruction_coefficients([ -1 ], k=self.k, ndiff=normal_derivative+1)

        src.append(set_coeffs(c, 0, stencil))

        if edge in [ 'right', 'left' ]:
            for s in stencil:
                src.append(set_vec('vec', [ phi[t,s] for t in stencil ]) 
                           + '; tmp(%+d) = dot_product(cof, vec)' % s)
        else:
            for s in stencil:
                src.append(set_vec('vec', [ phi[s,t] for t in stencil ]) 
                           + '; tmp(%+d) = dot_product(cof, vec)' % s)

        c = reconstruction_coefficients(self.quad_pts, k=self.k)

        for i in range(len(self.quad_pts)):
            src.append(set_coeffs(c, i, stencil))
            src.append('{unknown}_edge({qpt}) = dot_product(cof, tmp)'.format(
                unknown=dest, qpt=i))

        self.fsrc.extend(src)


    def reconstruct_interior(self, dest, unknown):
        """Generate code to reconstruct *unknown* in the interior."""

        if unknown in self.rcache:
            return
        else:
            self.rcache[unknown] = True

        src = []
        phi = make_unknown(unknown)
        stencil = self.stencil

        c = reconstruction_coefficients(self.quad_pts, k=self.k)

        for i in range(len(self.quad_pts)):
            src.append(set_coeffs(c, i, stencil))
            for s in stencil:
                src.append(set_vec('vec', [ phi[t,s] for t in stencil ]))
                src.append('tmp(%+d,%d) = dot_product(cof, vec)' % (s,i))

        for i in range(len(self.quad_pts)):
            for j in range(len(self.quad_pts)):
                src.append(set_coeffs(c, j, stencil))
                src.append('{unknown}_interior({i},{j}) = dot_product(cof, tmp(:,{i}))'.format(
                        unknown=dest, i=i, j=j))

        self.ssrc.extend(src)


    def add_edge_avg(self, dest, op, edge):

        self.rcache = {}

        # cycle through atoms in op and reconstruct edge values
        for name, unknown, ndiff in name_and_unknown(op):
            self.reconstruct_edge(name, unknown, edge, ndiff)

        # convert operator to code form
        self.fsrc.append('qvls = ' + to_string_op(op, '_edge'))
        
        if edge in [ 'right', 'top' ]:
            self.fsrc.append('{dest} = {dest} + 0.5d0*dx*dot_product(qwts, qvls)'.format(dest=dest))
        else:
            self.fsrc.append('{dest} = {dest} - 0.5d0*dx*dot_product(qwts, qvls)'.format(dest=dest))


    def add_vol_avg(self, dest, op):

        self.rcache = {}

        # cycle through atoms in op and reconstruct interior values
        for name, unknown, ndiff in name_and_unknown(op):
            self.reconstruct_interior(name, unknown)

        # convert operator to code form
        self.ssrc.append('qvls2 = ' + to_string_op(op, '_interior'))

        for i in range(len(self.quad_pts)):
            self.ssrc.append('qvls1({i}) = 0.5d0*dx*dot_product(qwts, qvls2(:, {i}))'.format(i=i))

        self.ssrc.append('{dest} = {dest} + 0.5d0*dx*dot_product(qwts, qvls1)'.format(dest=dest))


    def gen_flux(self, op, dest='f(i, j)'):
        """Generate code to compute the flux integral of *op*."""

        try:
            xop = op[0]
            yop = op[1]
        except:
            raise ValueError("Expected a vector function with two components.")

        self.set_quads()

        self.add_edge_avg(dest, op[0], 'right')
        self.add_edge_avg(dest, op[0], 'left')
        self.add_edge_avg(dest, op[1], 'top')
        self.add_edge_avg(dest, op[1], 'bottom')


    def flux(self):
        return '\n'.join(self.fsrc)


    def gen_source(self, op, dest='f(i, j)'):
        """Generate code to compute volume integral of *op*."""

        self.set_quads()

        self.add_vol_avg(dest, op)


    def source(self):
        return '\n'.join(self.ssrc)



