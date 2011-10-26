
import numpy as np
from pyfboxlib import fboxlib
import base

from fab import fab

class multifab(base.BLObject):
  """MultiFAB."""

  def create(self, layout, components=1, ghost_cells=0):
    """Create a multifab from a layout."""

    create   = getattr(fboxlib, 'create_' + self.__class__.__name__ + '_from_layout')
    get_info = getattr(fboxlib, 'get_' + self.__class__.__name__ + '_info')

    self.oid = create(layout.oid, components, ghost_cells)

    if self.oid:
      self.dim, self.nboxes, self.nc, self.ng = get_info(self.oid)

  def fab(self, box):
    return fab(self, box)


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









