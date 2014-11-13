"""PyBoxLib layout class."""

import numpy as np
import fcboxlib

class layout():
  """BoxLib layout."""

  def __init__(self, boxarray=None, cptr=None):
    if cptr:
      self.cptr = cptr
    else:
      pmask = np.ones(boxarray.dim, np.int32)
      self.cptr = fcboxlib.layout_create_from_boxarray(boxarray.cptr, pmask)

  @property
  def nboxes(self):
    """Return number of boxes."""
    return fcboxlib.layout_nboxes(self.cptr)

  def echo(self):
    fcboxlib.layout_print(self.cptr)

  # @property
  # def local_boxes(self):
  #   """Return list of local boxes."""
  #   boxes = []
  #   for k in range(self.nboxes):
  #     if self.local(k+1):
  #       boxes.append(k+1)
  #   return boxes


  # def local(self, box):
  #   """Return True if box is local."""
  #   local = c_int()
  #   bl.pybl_layout_local(self.cptr, box, byref(local))
  #   return local.value != 0


  # def get_box(self, i):
  #   bx = box()
  #   bl.pybl_layout_get_box(self.cptr, byref(bx), i)
  #   return bx


  # def from_regrid(self, lmfab):
  #   """Create a layout from re-gridding a tagged multifab."""

  #   ba = boxarray()
  #   ba.from_regrid(lmfab)

  #   self.create(boxarray=ba)
