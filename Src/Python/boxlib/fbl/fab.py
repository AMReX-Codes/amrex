"""PyBoxLib fab class."""

__all__ = [ 'fab' ]

from ctypes import *
from base import BLObject
from pybl import bl

# XXX: it would be nice to put libpycboxlib into libpyfboxlib
import libpycboxlib as cbl

class fab(BLObject):
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

  def __init__(self, mfab, nbox, level=None):

    self.mfab   = mfab
    self.nbox   = nbox
    self.bx_lo  = (3*c_int)()
    self.bx_hi  = (3*c_int)()
    self.pbx_lo = (3*c_int)()
    self.pbx_hi = (3*c_int)()
    self.ibx_lo = (3*c_int)()
    self.ibx_hi = (3*c_int)()

    self.c_int_attrs = [ 'dim', 'nc' ]
    self.init_c_int_attrs()

    assert 1 <= nbox <= mfab.nboxes

    mftype    = self.mfab.__class__.__name__
    get_info  = getattr(bl, 'pybl_get_' + mftype + '_fab_info')
    get_array = getattr(cbl, mftype + '_as_numpy')

    if level:
      get_info(self.mfab.cptr, level, nbox, byref(self._dim), byref(self._nc),
               self.pbx_lo, self.pbx_hi)
    else:
      get_info(self.mfab.cptr, nbox, byref(self._dim), byref(self._nc),
               self.bx_lo, self.bx_hi,
               self.pbx_lo, self.pbx_hi,
               self.ibx_lo, self.ibx_hi)

    ints = lambda arr: [ int(x) for x in arr[:self.dim] ]
    lohi = lambda lo, hi: (ints(lo), ints(hi))

    self.bx  = lohi(self.bx_lo,  self.bx_hi)
    self.ibx = lohi(self.ibx_lo, self.ibx_hi)
    self.pbx = lohi(self.pbx_lo, self.pbx_hi)

    cptr = addressof(self.mfab.cptr)
    if level:
      self.array = get_array(cptr, level, nbox).squeeze()
    else:
      self.array = get_array(cptr, nbox).squeeze()


  @property
  def shape(self):
    return self.array.shape

  @property
  def size(self):
    return self.array.size


  def bxrange(self, dim):
    """Return an iterator to cycle over the global indicies for the
    given dimension.

    Essentially this returns ``range(lo(dim), hi(dim)+1)``.
    """

    return range(self.bx[0][dim-1], self.bx[1][dim-1]+1)


  def pbxrange(self, dim):
    """Return an iterator to cycle over the global indicies for the
    given dimension.

    Essentially this returns ``range(lo(dim), hi(dim)+1)``.
    """

    return range(self.pbx[0][dim-1], self.pbx[1][dim-1]+1)


  def __getitem__(self, key):
    lbound = list(self.pbx[0])
    if len(self.array.shape) > self.dim:
      lbound.append(0)

    key = adjust_indexes(lbound, key)

    return self.array[key]


  def __setitem__(self, key, value):
    lbound = list(self.pbx[0])
    if len(self.array.shape) > self.dim:
      lbound.append(0)

    key = adjust_indexes(lbound, key)

    self.array[key] = value


# adapted from petsc4py
def adjust_indexes(lbounds, index):
   if not isinstance(index, tuple):
       return adjust_index(lounds[0], index)

   index = list(index)
   for i, start in enumerate(lbounds):
       index[i] = adjust_index(start, index[i])
   index = tuple(index)

   return index


# adapted from petsc4py
def adjust_index(lbound, index):

  if index is None:
      return index

  if index is Ellipsis:
      return index

  if isinstance(index, slice):
      start = index.start
      stop  = index.stop
      step  = index.step
      if start is not None: start -= lbound
      if stop  is not None: stop  -= lbound
      return slice(start, stop, step)

  try:
      return index - lbound

  except TypeError:
      return index
