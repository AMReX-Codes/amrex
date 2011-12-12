
import pyfboxlib
from pyfboxlib import fboxlib


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

    assert 1 <= nbox <= mfab.nboxes

    get_info = getattr(fboxlib, 'get_' + self.mfab.__class__.__name__ + '_fab_info')
    get_array = getattr(pyfboxlib, self.mfab.__class__.__name__ + '_array')

    (self.dim, self.nc, self.bx_lo, self.bx_hi, self.pbx_lo, self.pbx_hi,
     self.ibx_lo, self.ibx_hi) = get_info(mfab.oid, nbox)

    intarray = lambda a: [ int(e) for e in a[:self.dim] ]
    lohi = lambda lo, hi: (intarray(lo), intarray(hi))

    self.bx  = lohi(self.bx_lo,  self.bx_hi)
    self.ibx = lohi(self.ibx_lo, self.ibx_hi)
    self.pbx = lohi(self.pbx_lo, self.pbx_hi)

    self.array = get_array(self.mfab.oid, nbox).squeeze()

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
