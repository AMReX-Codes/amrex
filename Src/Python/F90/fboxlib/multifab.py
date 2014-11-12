"""PyBoxLib multifab class."""

from ctypes import *

import numpy as np
import base

from pybl import bl
from fab import fab

class multifab(base.BLObject):
  """MultiFAB."""

  def __init__(self, *args, **kwargs):

    super(multifab, self).__init__(*args, **kwargs)

    self.c_int_attrs = [ 'dim', 'nboxes', 'nc', 'ng' ]
    self.init_c_int_attrs()

  def __del__(self):
    bl.pybl_multifab_destroy(self.cptr)

  def create(self, layout, components=1, ghost_cells=0):
    """Create a multifab from a layout."""

    # XXX: nodal?

    mftype = self.__class__.__name__
    create = getattr(bl, 'pybl_create_' + mftype + '_from_layout')
    create(layout.cptr, components, ghost_cells, byref(self.cptr))

    self.get_info()


  def get_info(self):

    mftype = self.__class__.__name__

    if self.associated:
      get_info = getattr(bl, 'pybl_get_' + mftype + '_info')
      get_info(self.cptr,
               byref(self._dim),
               byref(self._nboxes),
               byref(self._nc),
               byref(self._ng))


  def create_from_bbox(self, mf, components=1, ghost_cells=0):
    """Creat a multifab from the bounding box of the existing mf multifab."""

    mftype = self.__class__.__name__

    create    = getattr(bl, 'pybl_create_' + mftype + '_from_bbox')
    self.cptr = create(mf.cptr, components, ghost_cells, byref(self.cptr))

    self.get_info()


  def copy(self, dest):
    """Copy self into *dest* (multifab)."""
    bl.pybl_multifab_copy(dest.cptr, self.cptr)


  def fill_boundary(self):
    """Fill ghost cells between processors."""
    bl.pybl_multifab_fill_boundary(self.cptr)


  def fab(self, i):
    """Return fab corresponding to box *i*."""
    return fab(self, i)


  def write(self, dirname, header):
    bl.pybl_multifab_write(self.cptr,
                           dirname, len(dirname),
                           header, len(header))


  def read(self, dirname, header):
    bl.pybl_multifab_read(dirname, len(dirname), header, len(header),
                          byref(self.cptr))
    self.get_info()


class lmultifab(multifab):
  """Logical MultiFAB."""

  pass
