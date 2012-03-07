"""PyBoxLib fab class."""

__all__ = [ 'fab' ]

from ctypes import *

import libpycboxlib as cbl
from pybl import bl


class fab(object):
  """FAB.

  Here we assume that a fab is attached to a multifab.

  DIM:    Dimension
  BX:     The Box in index space for which this FAB is defined
  IBX:    The index range of the valid data for the FAB
  PBX:    The physical box for the FAB
  NC:     Number of components
  NG:     Number of ghost cells

  When a FAB is created IBX = BX, unless it is nodal, in which case
  IBX = grow(BX, FACE=hi, 1 (in the nodal directions).
  PBX = grow(IBX, NG)

  For parallel systems the BX, IBX, PBX, etc are all defined, but the
  underlying pointer will not be allocated.

  All FABS are 'Four' dimensional, conventially, (NX,NY,NZ,NC) in size.
  NY = 1, NZ = 1, when DIM =1, NZ = 1, when DIM = 2.

  """

  def __init__(self, mfab, nbox):

    self.mfab = mfab
    self.nbox = nbox
    self.dim = c_int(0)
    self.nc = c_int(0)
    self.bx_lo = (3*c_int)()
    self.bx_hi = (3*c_int)()
    self.pbx_lo = (3*c_int)()
    self.pbx_hi = (3*c_int)()
    self.ibx_lo = (3*c_int)()
    self.ibx_hi = (3*c_int)()

    assert 1 <= nbox <= mfab.nboxes

    mftype = self.mfab.__class__.__name__
    get_info = getattr(bl, 'pybl_get_' + mftype + '_fab_info')
    get_array = getattr(cbl, mftype + '_array')

    get_info(self.mfab.cptr, nbox, byref(self.dim), byref(self.nc), 
             self.bx_lo, self.bx_hi, 
             self.pbx_lo, self.pbx_hi, 
             self.ibx_lo, self.ibx_hi)

    ints = lambda arr: [ int(x) for x in arr[:self.dim.value] ]
    lohi = lambda lo, hi: (ints(lo), ints(hi))

    self.bx  = lohi(self.bx_lo,  self.bx_hi)
    self.ibx = lohi(self.ibx_lo, self.ibx_hi)
    self.pbx = lohi(self.pbx_lo, self.pbx_hi)

    cptr = addressof(self.mfab.cptr)

    self.array = get_array(cptr, nbox).squeeze()


  @property
  def shape(self):
    return self.array.shape

  @property
  def size(self):
    return self.array.size

  def __getitem__(self, key):
    # XXX: switch to global indexing
    return self.array[key]

  def __setitem__(self, key, value):
    # XXX: switch to global indexing
    self.array[key] = value
