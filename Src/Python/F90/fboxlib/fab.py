"""PyBoxLib fab class."""

__all__ = [ 'fab' ]

import fboxlib.fcboxlib as fcboxlib
import numpy as np

class fab():
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

  def __init__(self, mfab, nbox, logical=False, squeeze=False):
    if logical:
      a = fcboxlib.lmultifab_as_numpy(mfab.cptr, nbox)
    else:
      a = fcboxlib.multifab_as_numpy(mfab.cptr, nbox)
    if squeeze:
      self.array = a.squeeze()
    else:
      shp = a.shape[:mfab.dim] + (a.shape[3],)
      self.array = a.reshape(shp)

  @property
  def shape(self):
    return self.array.shape

  @property
  def size(self):
    return self.array.size
