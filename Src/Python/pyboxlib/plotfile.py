"""PyBoxLib plotfile class."""

from ctypes import *

import base

from pybl import bl
from fab import fab


class plotfile(base.BLObject):
  """PlotFile."""


  def __init__(self, *args, **kwargs):

    super(plotfile, self).__init__(*args, **kwargs)

    self.c_int_attrs = [ 'dim', 'nvars', 'flevel' ]
    self.init_c_int_attrs()


  def create(self, dname):
    """Create a plotfile given a name."""

    create = getattr(bl, 'pybl_create_plotfile')
    create(dname, len(dname), byref(self.cptr))

    self.get_info()


  def get_info(self):

    # XXX: dx

    self.variables = {}

    if self.associated:
      get_info = getattr(bl, 'pybl_get_plotfile_info')
      get_info(self.cptr,
               byref(self._dim),
               byref(self._nvars),
               byref(self._flevel))

      get_name = getattr(bl, 'pybl_get_plotfile_name')
      for i in range(self.nvars):
        buf = create_string_buffer(20)
        get_name(self.cptr, i+1, 20, buf)
        self.variables[buf.value.strip()] = i+1


  def nboxes(self, level):
    nboxes = c_int()
    bl.pybl_get_plotfile_grid_info(self.cptr, level, byref(nboxes))
    return nboxes.value


  def bind(self, i, j, c):
    bl.pybl_plotfile_bind(self.cptr, i, j, c)


  def unbind(self, i, j):
    bl.pybl_plotfile_unbind(self.cptr, i, j)


  def fab(self, i, j):
    return fab(self, j, level=i)


def compare(dname1, dname2):
  """Compare two plotfiles."""

  p1 = plotfile()
  p1.create(dname1)

  p2 = plotfile()
  p2.create(dname2)

  if p1.dim != p2.dim:
    raise ValueError("Number of dimensions don't match: %s, %s" % (dname1, dname2))

  if p1.flevel != p2.flevel:
    raise ValueError("Number of levels don't match: %s, %s" % (dname1, dname2))

  errors = {}
  
  for variable in sorted(p1.variables):

    errors[variable] = {}

    c1 = p1.variables[variable]
    try:
      c2 = p2.variables[variable]
    except KeyError:
      print "WARNING: variable %s not found in plotfile %s, skipping..." % (variable, dname2)
      continue

    aerror = 0.0
    rerror = 0.0

    for i in range(1, p1.flevel+1):
      for j in range(1, p1.nboxes(i)+1):
        p1.bind(i, j, c1)
        p2.bind(i, j, c2)

        a1 = p1.fab(i, j)
        a2 = p2.fab(i, j)

        aerror += abs(a1.array - a2.array).max()
        rerror += aerror / abs(a1.array).max()

        p1.unbind(i, j)
        p2.unbind(i, j)

    errors[variable] = (aerror, rerror)

  return errors, (dname1, dname2)

def print_compare(errors, plotfiles):

  dname1, dname2 = plotfiles

  print ''
  print '=== error between %s and %s ===' % (dname1, dname2)
  print ''
  print '%20s %20s %20s' % ('', 'abs err', 'rel err')

  for variable in sorted(errors):
    aerror, rerror = errors[variable]

    print '%20s %20e %20e' % (variable, aerror, rerror)

  print ''
