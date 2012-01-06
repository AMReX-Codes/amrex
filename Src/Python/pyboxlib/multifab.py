
import numpy as np
from pyfboxlib import fboxlib
import base

from fab import fab

class multifab(base.BLObject):
  """MultiFAB."""

  def create(self, layout, components=1, ghost_cells=0, interleave=False):
    """Create a multifab from a layout."""

    create   = getattr(fboxlib, 'create_' + self.__class__.__name__ + '_from_layout')
    get_info = getattr(fboxlib, 'get_' + self.__class__.__name__ + '_info')

    self.oid = create(layout.oid, components, ghost_cells, interleave)
    self.interleaved = interleave

    if self.associated:
      self.dim, self.nboxes, self.nc, self.ng = get_info(self.oid)

  def create_from_bbox(self, mf, components=1, ghost_cells=0, interleave=False):
    """Creat a multifab from the bounding box of the existing mf multifab."""

    create   = getattr(fboxlib, 'create_' + self.__class__.__name__ + '_from_bbox')
    get_info = getattr(fboxlib, 'get_' + self.__class__.__name__ + '_info')

    self.oid = create(mf.oid, components, ghost_cells, interleave)
    self.interleaved = interleave

    if self.associated:
      self.dim, self.nboxes, self.nc, self.ng = get_info(self.oid)

  def copy(self, dmfab):
    """Copy self into dmfab."""
    fboxlib.pybl_multifab_copy(dmfab.oid, self.oid)

  def fab(self, box):
    return fab(self, box)

  def write(self, dirname, header):
    fboxlib.pybl_multifab_write(self.oid, dirname, header)

  def read(self, dirname, header):
    self.oid = fboxlib.pybl_multifab_read(dirname, header)

class lmultifab(multifab):
  """Logical MultiFAB."""

  pass

  # def create(self, layout):
  #   """Create a logical (boolean) multifab from a layout."""

  #   self.oid = fboxlib.create_lmultifab_from_layout(layout.oid)

  #   if self.oid:
  #     self.dim, self.nboxes = fboxlib.get_lmultifab_info(self.oid)

  # # # XXX: add a more fancy get/set
  # # def array(self, box):
  # #   return lmultifab_array(self.oid, box)









