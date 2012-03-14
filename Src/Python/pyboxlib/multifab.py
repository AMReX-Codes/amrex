"""PyBoxLib multifab class."""

from ctypes import *

import numpy as np
import base

from pybl import bl
from fab import fab

class multifab(base.BLObject):
  """MultiFAB."""

  def create(self, layout, components=1, ghost_cells=0, interleave=False):
    """Create a multifab from a layout."""

    self.dim = c_int(0)
    self.nboxes = c_int(0)
    self.nc = c_int(0)
    self.ng = c_int(0)
    self.interleaved = c_int(0)

    mftype = self.__class__.__name__
    create = getattr(bl, 'pybl_create_' + mftype + '_from_layout')
    create(layout.cptr, components, ghost_cells, interleave, byref(self.cptr))

    if self.associated:
      get_info = getattr(bl, 'pybl_get_' + mftype + '_info')
      get_info(self.cptr,
               byref(self.dim), 
               byref(self.nboxes), 
               byref(self.nc), 
               byref(self.ng))


  def create_from_bbox(self, mf, components=1, ghost_cells=0, interleave=False):
    """Creat a multifab from the bounding box of the existing mf multifab."""

    mftype = self.__class__.__name__ 

    create   = getattr(bl, 'create_' + mftype + '_from_bbox')
    get_info = getattr(bl, 'get_' + mftype + '_info')

    self.cptr = create(mf.cptr, components, ghost_cells, interleave)
    self.interleaved = interleave

    if self.associated:
      self.dim, self.nboxes, self.nc, self.ng = get_info(self.cptr)


  def copy(self, dest):
    """Copy self into *dest* (multifab)."""

    bl.pybl_multifab_copy(dest.cptr, self.cptr)


  def fill_boundary(self):
    """Fill ghost cells between processors."""

    bl.pybl_multifab_fill_boundary(self.cptr)


  def fab(self, box):
    """Return fab corresponding to box *box*."""

    return fab(self, box)


  def write(self, dirname, header):
    bl.pybl_multifab_write(self.cptr, dirname, header)


  def read(self, dirname, header):
    self.cptr = bl.pybl_multifab_read(dirname, header)


class lmultifab(multifab):
  """Logical MultiFAB."""

  pass

  # def create(self, layout):
  #   """Create a logical (boolean) multifab from a layout."""

  #   self.cptr = pybl.create_lmultifab_from_layout(layout.cptr)

  #   if self.cptr:
  #     self.dim, self.nboxes = pybl.get_lmultifab_info(self.cptr)

  # # # XXX: add a more fancy get/set
  # # def array(self, box):
  # #   return lmultifab_array(self.cptr, box)









